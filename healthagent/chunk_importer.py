"""Chunked DNA file importer using threading, Redis, and MongoDB.

Pipeline:
  1. SystemProbe  — measure available CPU / RAM and compute optimal chunk config
  2. DNAChunkImporter.import_file() — split raw file into chunks, dispatch to
     a ThreadPoolExecutor, store parsed SNPs in per-chromosome MongoDB
     collections, track progress in Redis.

Redis key schema:
  dna:job:<job_id>:status        STRING  "queued" | "running" | "done" | "error"
  dna:job:<job_id>:total_chunks  STRING  total number of chunks
  dna:job:<job_id>:done_chunks   STRING  atomic counter incremented per chunk
  dna:job:<job_id>:snp_count     STRING  total SNPs inserted
  dna:job:<job_id>:chunks        LIST    chunk metadata JSON entries

MongoDB collection schema (one collection per chromosome):
  db: healthagent
  collection: chr_<chromosome>   e.g. chr_1, chr_X, chr_MT
  document: { rsid, chromosome, position, genotype, job_id, chunk_index }
"""

import json
import os
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import psutil
import redis
from pymongo import MongoClient, ASCENDING
from pymongo.collection import Collection

from healthagent.dna_importer import DNAFormat, SNP, import_dna_string


# ---------------------------------------------------------------------------
# System probe
# ---------------------------------------------------------------------------

@dataclass
class SystemResources:
    cpu_cores: int
    cpu_percent_free: float
    ram_total_mb: int
    ram_available_mb: int
    recommended_threads: int
    recommended_chunk_lines: int

    def report(self) -> str:
        return (
            f"CPU cores: {self.cpu_cores}  "
            f"CPU free: {self.cpu_percent_free:.1f}%  "
            f"RAM total: {self.ram_total_mb} MB  "
            f"RAM available: {self.ram_available_mb} MB  "
            f"=> threads: {self.recommended_threads}  "
            f"chunk size: {self.recommended_chunk_lines:,} lines"
        )


def probe_system(
    min_threads: int = 2,
    max_threads: int = 32,
    min_chunk_lines: int = 5_000,
    max_chunk_lines: int = 500_000,
) -> SystemResources:
    """Measure live CPU and RAM and compute safe threading / chunk parameters.

    Rules:
    - Use at most 75 % of available CPU cores (leave headroom for OS).
    - Each worker thread is allocated ~64 MB of RAM budget; scale down if RAM
      is tight (< 512 MB available → fall back to 2 threads, small chunks).
    - Chunk size scales linearly between min/max based on available RAM.
    """
    cpu_cores = psutil.cpu_count(logical=True) or 1
    cpu_pct_used = psutil.cpu_percent(interval=0.5)
    vm = psutil.virtual_memory()
    ram_total_mb = vm.total // (1024 ** 2)
    ram_avail_mb = vm.available // (1024 ** 2)

    # Thread budget: 75% of cores, capped by RAM (64 MB per thread)
    threads_by_cpu = max(min_threads, int(cpu_cores * 0.75))
    threads_by_ram = max(min_threads, ram_avail_mb // 64)
    threads = min(threads_by_cpu, threads_by_ram, max_threads)

    # Chunk size: scale linearly with available RAM
    ram_fraction = min(1.0, ram_avail_mb / 4096)  # saturates at 4 GB
    chunk_lines = int(
        min_chunk_lines + ram_fraction * (max_chunk_lines - min_chunk_lines)
    )

    return SystemResources(
        cpu_cores=cpu_cores,
        cpu_percent_free=100.0 - cpu_pct_used,
        ram_total_mb=ram_total_mb,
        ram_available_mb=ram_avail_mb,
        recommended_threads=threads,
        recommended_chunk_lines=chunk_lines,
    )


# ---------------------------------------------------------------------------
# Chunked importer
# ---------------------------------------------------------------------------

@dataclass
class ChunkResult:
    chunk_index: int
    snp_count: int
    chromosomes: list[str]
    elapsed_sec: float
    error: Optional[str] = None


class DNAChunkImporter:
    """Import a raw DNA file in parallel chunks into MongoDB, tracked via Redis.

    Args:
        mongo_uri:   MongoDB connection string (default: localhost).
        redis_host:  Redis host (default: localhost).
        redis_port:  Redis port (default: 6379).
        db_name:     MongoDB database name (default: "healthagent").
    """

    def __init__(
        self,
        mongo_uri: str = "mongodb://localhost:27017/",
        redis_host: str = "localhost",
        redis_port: int = 6379,
        db_name: str = "healthagent",
    ):
        self.mongo_uri = mongo_uri
        self.redis_host = redis_host
        self.redis_port = redis_port
        self.db_name = db_name
        self._lock = threading.Lock()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def import_file(
        self,
        path: str | Path,
        job_id: Optional[str] = None,
        resources: Optional[SystemResources] = None,
    ) -> dict:
        """Import a DNA file using adaptive threading and chunking.

        Args:
            path:      Path to raw DNA file.
            job_id:    Optional stable job identifier (UUID generated if omitted).
            resources: Pre-computed SystemResources; probed live if omitted.

        Returns:
            Summary dict with job_id, snp_count, chunk_count, collections, elapsed.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"DNA file not found: {path}")

        job_id = job_id or str(uuid.uuid4())
        resources = resources or probe_system()

        print(f"[job {job_id}] {resources.report()}")

        raw = path.read_text(encoding="utf-8", errors="replace")
        return self._run(raw, job_id, resources, filename=path.name)

    def import_string(
        self,
        raw: str,
        job_id: Optional[str] = None,
        resources: Optional[SystemResources] = None,
        filename: str = "",
    ) -> dict:
        """Import DNA data from a raw string."""
        job_id = job_id or str(uuid.uuid4())
        resources = resources or probe_system()
        return self._run(raw, job_id, resources, filename=filename)

    # ------------------------------------------------------------------
    # Internal pipeline
    # ------------------------------------------------------------------

    def _run(self, raw: str, job_id: str, res: SystemResources, filename: str) -> dict:
        t_start = time.monotonic()

        # Split into data lines (skip blank + comment lines)
        all_lines = raw.splitlines()
        header_lines = [l for l in all_lines if l.startswith("#")]
        data_lines = [l for l in all_lines if l.strip() and not l.startswith("#")]

        # Detect format once using the full header, reconstruct per-chunk preamble
        from healthagent.dna_importer import detect_format
        first_data = data_lines[0] if data_lines else ""
        fmt = detect_format(header_lines, first_data)
        preamble = "\n".join(header_lines) + "\n"  # prepended to every chunk

        chunks = self._split_chunks(data_lines, res.recommended_chunk_lines)
        total_chunks = len(chunks)

        # Redis: register job
        r = self._redis()
        r.set(f"dna:job:{job_id}:status", "running")
        r.set(f"dna:job:{job_id}:total_chunks", total_chunks)
        r.set(f"dna:job:{job_id}:done_chunks", 0)
        r.set(f"dna:job:{job_id}:snp_count", 0)
        r.set(f"dna:job:{job_id}:filename", filename)
        r.set(f"dna:job:{job_id}:format", fmt.value)

        print(
            f"[job {job_id}] {len(data_lines):,} data lines → "
            f"{total_chunks} chunks × ~{res.recommended_chunk_lines:,} lines  "
            f"using {res.recommended_threads} threads"
        )

        results: list[ChunkResult] = []
        collections_used: set[str] = set()

        with ThreadPoolExecutor(max_workers=res.recommended_threads) as pool:
            futures = {
                pool.submit(
                    self._process_chunk,
                    preamble + "\n".join(chunk),
                    idx,
                    job_id,
                ): idx
                for idx, chunk in enumerate(chunks)
            }
            for future in as_completed(futures):
                result: ChunkResult = future.result()
                results.append(result)
                collections_used.update(result.chromosomes)

                # Update Redis progress atomically
                r.incr(f"dna:job:{job_id}:done_chunks")
                r.incrby(f"dna:job:{job_id}:snp_count", result.snp_count)
                r.lpush(
                    f"dna:job:{job_id}:chunks",
                    json.dumps(
                        {
                            "chunk_index": result.chunk_index,
                            "snp_count": result.snp_count,
                            "chromosomes": result.chromosomes,
                            "elapsed_sec": round(result.elapsed_sec, 3),
                            "error": result.error,
                        }
                    ),
                )

                done = int(r.get(f"dna:job:{job_id}:done_chunks"))
                pct = done / total_chunks * 100
                print(
                    f"[job {job_id}] chunk {result.chunk_index:>4}/{total_chunks}  "
                    f"{pct:5.1f}%  SNPs={result.snp_count}  "
                    f"chrs={result.chromosomes}  {result.elapsed_sec:.2f}s"
                )

        errors = [r for r in results if r.error]
        r.set(f"dna:job:{job_id}:status", "done" if not errors else "partial")

        total_snps = sum(r.snp_count for r in results)
        elapsed = time.monotonic() - t_start

        summary = {
            "job_id": job_id,
            "status": "done" if not errors else "partial",
            "format": fmt.value,
            "filename": filename,
            "chunk_count": total_chunks,
            "snp_count": total_snps,
            "collections": sorted(f"chr_{c}" for c in collections_used),
            "threads_used": res.recommended_threads,
            "chunk_size_lines": res.recommended_chunk_lines,
            "elapsed_sec": round(elapsed, 3),
            "errors": [{"chunk": e.chunk_index, "msg": e.error} for e in errors],
        }
        print(f"[job {job_id}] DONE  {total_snps:,} SNPs in {elapsed:.2f}s")
        return summary

    # ------------------------------------------------------------------
    # Per-chunk worker (runs in thread)
    # ------------------------------------------------------------------

    def _process_chunk(self, chunk_raw: str, chunk_index: int, job_id: str) -> ChunkResult:
        t0 = time.monotonic()
        try:
            profile = import_dna_string(chunk_raw)
        except Exception as exc:
            return ChunkResult(chunk_index=chunk_index, snp_count=0, chromosomes=[], elapsed_sec=time.monotonic() - t0, error=str(exc))

        if not profile.snps:
            return ChunkResult(chunk_index=chunk_index, snp_count=0, chromosomes=[], elapsed_sec=time.monotonic() - t0)

        # Group SNPs by chromosome → separate MongoDB collections
        by_chrom: dict[str, list[dict]] = {}
        for snp in profile.snps:
            doc = {
                "rsid": snp.rsid,
                "chromosome": snp.chromosome,
                "position": snp.position,
                "genotype": snp.genotype,
                "job_id": job_id,
                "chunk_index": chunk_index,
            }
            by_chrom.setdefault(snp.chromosome, []).append(doc)

        # Each thread gets its own MongoDB client (not thread-safe to share)
        client = MongoClient(self.mongo_uri, serverSelectionTimeoutMS=5000)
        db = client[self.db_name]

        for chrom, docs in by_chrom.items():
            col: Collection = db[f"chr_{chrom}"]
            col.insert_many(docs, ordered=False)
            # Ensure index on rsid + position for fast lookups
            col.create_index([("rsid", ASCENDING)], background=True)
            col.create_index([("position", ASCENDING)], background=True)

        client.close()

        return ChunkResult(
            chunk_index=chunk_index,
            snp_count=len(profile.snps),
            chromosomes=list(by_chrom.keys()),
            elapsed_sec=time.monotonic() - t0,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _split_chunks(lines: list[str], chunk_size: int) -> list[list[str]]:
        """Split a list of lines into fixed-size chunks."""
        return [lines[i : i + chunk_size] for i in range(0, len(lines), chunk_size)]

    def _redis(self) -> redis.Redis:
        return redis.Redis(
            host=self.redis_host,
            port=self.redis_port,
            decode_responses=True,
        )

    # ------------------------------------------------------------------
    # Job status helpers
    # ------------------------------------------------------------------

    def job_status(self, job_id: str) -> dict:
        """Return current Redis-tracked status for a job."""
        r = self._redis()
        total = int(r.get(f"dna:job:{job_id}:total_chunks") or 0)
        done = int(r.get(f"dna:job:{job_id}:done_chunks") or 0)
        snps = int(r.get(f"dna:job:{job_id}:snp_count") or 0)
        return {
            "job_id": job_id,
            "status": r.get(f"dna:job:{job_id}:status"),
            "filename": r.get(f"dna:job:{job_id}:filename"),
            "format": r.get(f"dna:job:{job_id}:format"),
            "total_chunks": total,
            "done_chunks": done,
            "progress_pct": round(done / total * 100, 1) if total else 0,
            "snp_count": snps,
        }

    def list_collections(self) -> list[str]:
        """Return all chromosome collections in the MongoDB database."""
        client = MongoClient(self.mongo_uri, serverSelectionTimeoutMS=5000)
        cols = sorted(client[self.db_name].list_collection_names())
        client.close()
        return cols
