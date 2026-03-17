"""Core DNA health agent module."""

from pathlib import Path

from healthagent.dna_importer import DNAProfile, import_dna_file, import_dna_string


class HealthAgent:
    """AI agent for DNA health and welfare analysis."""

    def __init__(self):
        self.name = "HealthAgent"
        self.profile: DNAProfile | None = None

    def load_file(self, path: str | Path) -> DNAProfile:
        """Load a DNA file from disk (auto-detects format)."""
        self.profile = import_dna_file(path)
        return self.profile

    def load_string(self, raw: str, filename: str = "") -> DNAProfile:
        """Load DNA data from a raw string (auto-detects format)."""
        self.profile = import_dna_string(raw, filename=filename)
        return self.profile

    def analyze(self, dna_data: str) -> dict:
        """Analyze DNA data and return health insights."""
        raise NotImplementedError("Analysis not yet implemented")
