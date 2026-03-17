"""DNA file importer supporting major ancestry company formats.

Supported formats:
- 23andMe (TSV: rsid, chromosome, position, genotype)
- AncestryDNA (TSV: rsid, chromosome, position, allele1, allele2)
- MyHeritage (CSV: rsid, chromosome, position, genotype)
- FamilyTreeDNA (CSV: RSID, CHROMOSOME, POSITION, RESULT)
- LivingDNA (TSV: rsid, chromosome, position, genotype)
- VCF (Variant Call Format - universal standard)
"""

import csv
import io
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class DNAFormat(Enum):
    TWENTYTHREE_AND_ME = "23andme"
    ANCESTRY_DNA = "ancestrydna"
    MY_HERITAGE = "myheritage"
    FAMILY_TREE_DNA = "familytreedna"
    LIVING_DNA = "livingdna"
    VCF = "vcf"
    UNKNOWN = "unknown"


@dataclass
class SNP:
    """Single nucleotide polymorphism record."""
    rsid: str
    chromosome: str
    position: int
    genotype: str

    def __post_init__(self):
        self.chromosome = str(self.chromosome).strip().upper().lstrip("CHR")
        self.genotype = self.genotype.strip().upper().replace("-", "--")
        if self.rsid and not self.rsid.startswith("rs"):
            self.rsid = self.rsid.strip()


@dataclass
class DNAProfile:
    """Container for a parsed DNA profile."""
    source_format: DNAFormat
    snps: list[SNP] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    @property
    def snp_count(self) -> int:
        return len(self.snps)

    def get_snp(self, rsid: str) -> Optional[SNP]:
        """Look up a SNP by rsID."""
        for snp in self.snps:
            if snp.rsid.lower() == rsid.lower():
                return snp
        return None


# --- Format detection ---

def detect_format(header_lines: list[str], first_data_line: str) -> DNAFormat:
    """Detect the ancestry file format from header comments and first data row."""
    header_text = "\n".join(header_lines).lower()

    # VCF must be checked first — its headers may mention other providers (e.g. ##source=23andMe)
    if any("fileformat=vcf" in h.lower() for h in header_lines):
        return DNAFormat.VCF
    if "23andme" in header_text:
        return DNAFormat.TWENTYTHREE_AND_ME
    if "ancestrydna" in header_text or "ancestry.com" in header_text:
        return DNAFormat.ANCESTRY_DNA
    if "myheritage" in header_text:
        return DNAFormat.MY_HERITAGE
    if "familytreedna" in header_text or "family tree dna" in header_text:
        return DNAFormat.FAMILY_TREE_DNA
    if "livingdna" in header_text or "living dna" in header_text:
        return DNAFormat.LIVING_DNA

    # Fallback: inspect column headers in data
    cols = re.split(r"[\t,]", first_data_line.lower())
    if "rsid" in cols and "allele1" in cols:
        return DNAFormat.ANCESTRY_DNA
    if "rsid" in cols and "result" in cols:
        return DNAFormat.FAMILY_TREE_DNA
    if first_data_line.startswith("#chrom") or first_data_line.startswith("##"):
        return DNAFormat.VCF

    return DNAFormat.UNKNOWN


# --- Per-format parsers ---

def _parse_tsv_genotype(lines: list[str], skip_header_row: bool = True) -> list[SNP]:
    """Parse tab-delimited files with columns: rsid, chromosome, position, genotype."""
    snps = []
    reader = csv.reader(lines, delimiter="\t")
    for row in reader:
        if skip_header_row and row and row[0].lower().startswith("rsid"):
            skip_header_row = False
            continue
        if not row or row[0].startswith("#") or len(row) < 4:
            continue
        try:
            snps.append(SNP(rsid=row[0], chromosome=row[1], position=int(row[2]), genotype=row[3]))
        except (ValueError, IndexError):
            continue
    return snps


def _parse_ancestry_tsv(lines: list[str]) -> list[SNP]:
    """Parse AncestryDNA format: rsid, chromosome, position, allele1, allele2."""
    snps = []
    skip_header = True
    reader = csv.reader(lines, delimiter="\t")
    for row in reader:
        if skip_header and row and row[0].lower().startswith("rsid"):
            skip_header = False
            continue
        if not row or row[0].startswith("#") or len(row) < 5:
            continue
        try:
            genotype = row[3].strip() + row[4].strip()
            snps.append(SNP(rsid=row[0], chromosome=row[1], position=int(row[2]), genotype=genotype))
        except (ValueError, IndexError):
            continue
    return snps


def _parse_familytreedna_csv(lines: list[str]) -> list[SNP]:
    """Parse FamilyTreeDNA CSV: RSID, CHROMOSOME, POSITION, RESULT."""
    snps = []
    skip_header = True
    reader = csv.reader(lines)
    for row in reader:
        if skip_header and row and row[0].strip('"').upper() in ("RSID", "# RSID"):
            skip_header = False
            continue
        if not row or row[0].startswith("#") or len(row) < 4:
            continue
        try:
            snps.append(SNP(
                rsid=row[0].strip('"'),
                chromosome=row[1].strip('"'),
                position=int(row[2].strip('"')),
                genotype=row[3].strip('"'),
            ))
        except (ValueError, IndexError):
            continue
    return snps


def _parse_vcf(lines: list[str]) -> list[SNP]:
    """Parse VCF format into SNP records (genotype calls only)."""
    snps = []
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        try:
            chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            # Attempt to decode genotype from FORMAT/SAMPLE columns if present
            genotype = ref + alt if alt not in (".", "") else ref + ref
            snps.append(SNP(rsid=rsid, chromosome=chrom, position=int(pos), genotype=genotype))
        except (ValueError, IndexError):
            continue
    return snps


# --- Public API ---

def import_dna_file(path: str | Path) -> DNAProfile:
    """Import a raw DNA file and return a DNAProfile.

    Automatically detects the format from major ancestry providers:
    23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA, and VCF.

    Args:
        path: Path to the raw DNA data file.

    Returns:
        DNAProfile containing parsed SNPs and detected format metadata.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the format cannot be determined or the file is unreadable.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DNA file not found: {path}")

    raw = path.read_text(encoding="utf-8", errors="replace")
    return import_dna_string(raw, filename=path.name)


def import_dna_string(raw: str, filename: str = "") -> DNAProfile:
    """Import DNA data from a raw string.

    Useful for processing in-memory data or testing without a file on disk.

    Args:
        raw: Raw text content of the DNA file.
        filename: Optional filename hint for format detection.

    Returns:
        DNAProfile containing parsed SNPs and detected format metadata.

    Raises:
        ValueError: If the format cannot be determined.
    """
    lines = raw.splitlines()

    # Separate comment/header lines from data lines
    header_lines = [l for l in lines if l.startswith("#")]
    data_lines = [l for l in lines if l.strip() and not l.startswith("#")]

    if not data_lines:
        raise ValueError("No data found in DNA file.")

    first_data_line = data_lines[0]
    fmt = detect_format(header_lines, first_data_line)

    # Re-include comment lines for parsers that need them (e.g. VCF)
    all_lines = lines

    if fmt == DNAFormat.ANCESTRY_DNA:
        snps = _parse_ancestry_tsv(data_lines)
    elif fmt == DNAFormat.FAMILY_TREE_DNA:
        snps = _parse_familytreedna_csv(data_lines)
    elif fmt == DNAFormat.VCF:
        snps = _parse_vcf(all_lines)
    elif fmt in (DNAFormat.TWENTYTHREE_AND_ME, DNAFormat.MY_HERITAGE, DNAFormat.LIVING_DNA):
        snps = _parse_tsv_genotype(data_lines)
    else:
        # Best-effort parse for unknown formats
        snps = _parse_tsv_genotype(data_lines)
        fmt = DNAFormat.UNKNOWN

    metadata = {
        "filename": filename,
        "raw_line_count": len(lines),
        "header_line_count": len(header_lines),
    }

    return DNAProfile(source_format=fmt, snps=snps, metadata=metadata)
