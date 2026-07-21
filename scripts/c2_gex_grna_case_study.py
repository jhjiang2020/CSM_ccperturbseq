#!/usr/bin/env python3
"""Build the C2 gRNA-to-GEX misalignment case study.

The ``extract`` command is intended for a bound UGER compute job.  It reads the
completed C2 controlled-alignment audit, streams the original R1 FASTQ once,
and streams the query-name-sorted production GEX BAM once through samtools.
Only compact TSV/JSON/FASTQ evidence is retained.

The ``render`` command consumes that compact bundle and produces repository
figures and a Markdown report.  It deliberately does not import pysam so the
report can be rendered in the lightweight local analysis environment.

The ``publish`` command validates the internal extraction bundle and emits the
deterministic scientific-minimum evidence that is safe to version with the
report. It removes direct identifiers and filesystem paths before publication.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import pandas as pd
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


TSO = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
SPLITPIPE_TSO = "AACGCAGAGTGAATGGG"
STEM_VARIANTS = (
    "GTTTGAGAGCTAAGCAGGA",
    "GTTTCAGAGCTATGCTGGA",
)
TRACR_VARIANTS = (
    "GACTGCTTAGCAAGTT",
    "AACAGCATAGCAAGTT",
)
IBAR_LENGTH = 7
R2_LINKER_2 = "TCCAACCACCTC"
R2_LINKER_2_RC = "GAGGTGGTTGGA"

COMPONENT_COLORS = {
    "TSO": "#4C78A8",
    "guide": "#D62728",
    "stem": "#F2A541",
    "ibar": "#9467BD",
    "tracr": "#59A14F",
    "bridge": "#B9BEC7",
    "r2_barcode_rc": "#8C564B",
    "r2_linker_rc": "#17BECF",
}

@dataclass(frozen=True)
class TargetLocus:
    gene: str
    gene_id: str
    star_contig: str
    fasta_contig: str
    start_1: int
    end_1: int
    expected_region: str
    expected_target_reads: int
    expected_removed_false_molecules: int


TARGETS: Tuple[TargetLocus, ...] = (
    TargetLocus(
        gene="KLF2",
        gene_id="ENSG00000127528",
        star_contig="hg38_chr19",
        fasta_contig="chr19",
        start_1=16_324_750,
        end_1=16_325_250,
        expected_region="intronic",
        expected_target_reads=567,
        expected_removed_false_molecules=484,
    ),
    TargetLocus(
        gene="NFKB1",
        gene_id="ENSG00000109320",
        star_contig="hg38_chr4",
        fasta_contig="chr4",
        start_1=102_501_300,
        end_1=102_502_200,
        expected_region="exonic",
        expected_target_reads=500,
        expected_removed_false_molecules=414,
    ),
)

EXPECTED_GLOBAL = {
    "input_reads": 10_559,
    "full_mapped_primary_reads": 0,
    "short_mapped_primary_reads": 8_113,
    "short_unique_reads": 6_965,
    "short_multimapping_primary_reads": 1_148,
    "short_false_gene_assigned_reads": 4_008,
    "short_guide_target_assigned_reads": 2_431,
}

PUBLISHED_BUNDLE_SCHEMA = 1
PUBLISHED_CONTROLLED_COLUMNS = (
    "truth_id",
    "guide_query",
    "target_genes",
    "production_filter_removed",
    "contig",
    "start_1",
    "end_1",
    "blocks",
    "cigar",
    "strand",
    "mapq",
    "nh",
    "nm",
    "assigned_gene_id",
    "assigned_gene_name",
    "exonic_fraction",
    "annotation_region",
    "target_concordant",
)
PUBLISHED_REPRESENTATIVE_COLUMNS = (
    "truth_id",
    "guide_query",
    "target_genes",
    "production_filter_removed",
    "contig",
    "start_1",
    "end_1",
    "blocks",
    "cigar",
    "strand",
    "mapq",
    "nh",
    "nm",
    "assigned_gene_id",
    "assigned_gene_name",
    "exonic_fraction",
    "annotation_region",
    "target_concordant",
    "raw_r1",
    "synthetic54_r1",
    "processed_short54_r1",
    "processed_full_r1",
    "raw_tso",
    "raw_guide",
    "raw_post_guide",
    "reference_aligned",
    "query_aligned",
    "match_symbols",
    "full_arm_mapped",
)
PUBLISHED_GEX_COLUMNS = (
    "gene",
    "gene_id",
    "read_id",
    "contig",
    "start_1",
    "end_1",
    "blocks",
    "cigar",
    "strand",
    "mapq",
    "nm",
    "nh",
    "gx",
    "gn",
    "query_length",
    "post_splitpipe_structure_truth",
)


def normalize_gene_name(value: object) -> str:
    """Strip repeated GTF serialization quotes from a gene name."""

    name = str(value or "").strip()
    while len(name) >= 2 and name[0] == name[-1] and name[0] in {'"', "'"}:
        name = name[1:-1].strip()
    return name


def as_bool(value: object) -> bool:
    """Interpret booleans consistently after TSV/CSV serialization."""

    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def parse_truth_qname(query_name: str) -> Dict[str, str]:
    parts = query_name.split("|")
    parsed = {"truth_id": parts[0]}
    for field in parts[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            parsed[key] = value
    return parsed


def parse_cigar_blocks(start_1: int, cigar: str) -> List[Tuple[int, int]]:
    """Return 1-based inclusive aligned reference blocks for a SAM CIGAR."""

    if not cigar or cigar == "*":
        return []
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not tokens or "".join(f"{n}{op}" for n, op in tokens) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    position = int(start_1)
    block_start = None
    blocks: List[Tuple[int, int]] = []
    for length_text, op in tokens:
        length = int(length_text)
        if op in "M=X":
            if block_start is None:
                block_start = position
            position += length
        elif op in "DN":
            if block_start is not None:
                blocks.append((block_start, position - 1))
                block_start = None
            position += length
        elif op in "ISHP":
            continue
    if block_start is not None:
        blocks.append((block_start, position - 1))
    return blocks


def serialize_blocks(blocks: Sequence[Tuple[int, int]]) -> str:
    return ";".join(f"{start}-{end}" for start, end in blocks)


def deserialize_blocks(value: object) -> List[Tuple[int, int]]:
    text = str(value or "")
    if not text or text == "nan":
        return []
    blocks = []
    for token in text.split(";"):
        start, end = token.split("-", 1)
        blocks.append((int(start), int(end)))
    return blocks


def intervals_overlap(
    blocks: Sequence[Tuple[int, int]], start_1: int, end_1: int
) -> bool:
    return any(start <= end_1 and end >= start_1 for start, end in blocks)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted((int(start), int(end)) for start, end in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def splitpipe_trim(sequence: str) -> Tuple[str, int]:
    index = sequence.find(SPLITPIPE_TSO, 0, 15 + len(SPLITPIPE_TSO))
    if index < 0:
        return sequence, -1
    return sequence[index + len(SPLITPIPE_TSO) :], index


def raw_component_segments(
    raw_sequence: str, guide_sequence: str
) -> List[Tuple[int, int, str]]:
    """Resolve CSM and opposite-end read-through components in a raw R1.

    Coordinates are zero-based half-open intervals. The exact Parse-v3 second
    R2 linker anchors the opposite-end region; the eight bases immediately on
    either side are the reverse-oriented BC1 and BC2 fields, respectively.
    """

    raw = str(raw_sequence).upper()
    guide = str(guide_sequence).upper()
    if not raw.startswith(TSO + guide):
        raise ValueError("Raw representative does not start with exact TSO and guide")

    guide_start = len(TSO)
    stem_start = guide_start + len(guide)
    variant = next(
        (index for index, stem in enumerate(STEM_VARIANTS) if raw.startswith(stem, stem_start)),
        None,
    )
    if variant is None:
        raise ValueError("Raw representative lacks an exact configured CSM stem")

    stem_end = stem_start + len(STEM_VARIANTS[variant])
    ibar_end = stem_end + IBAR_LENGTH
    tracr_end = ibar_end + len(TRACR_VARIANTS[variant])
    if tracr_end > len(raw):
        raise ValueError("Raw representative ends before the partial tracr region")

    segments: List[Tuple[int, int, str]] = [
        (0, guide_start, "TSO"),
        (guide_start, stem_start, "guide"),
        (stem_start, stem_end, "stem"),
        (stem_end, ibar_end, "ibar"),
        (ibar_end, tracr_end, "tracr"),
    ]
    linker_start = raw.find(R2_LINKER_2_RC, tracr_end)
    if linker_start < 0:
        segments.append((tracr_end, len(raw), "bridge"))
        return segments

    bc1_start = linker_start - 8
    if bc1_start < tracr_end:
        raise ValueError("R2 linker RC is too close to the partial tracr")
    if bc1_start > tracr_end:
        segments.append((tracr_end, bc1_start, "bridge"))
    segments.extend(
        [
            (bc1_start, linker_start, "r2_barcode_rc"),
            (linker_start, linker_start + len(R2_LINKER_2_RC), "r2_linker_rc"),
        ]
    )
    bc2_start = linker_start + len(R2_LINKER_2_RC)
    bc2_end = min(len(raw), bc2_start + 8)
    if bc2_end > bc2_start:
        segments.append((bc2_start, bc2_end, "r2_barcode_rc"))
    if bc2_end < len(raw):
        segments.append((bc2_end, len(raw), "r2_linker_rc"))
    return segments


def is_exact_construct(
    raw_sequence: str, guide_sequence: str, stems: Sequence[str]
) -> bool:
    if len(raw_sequence) != 150 or not raw_sequence.startswith(TSO):
        return False
    post_tso = raw_sequence[len(TSO) :]
    if not post_tso.startswith(guide_sequence):
        return False
    post_guide = post_tso[len(guide_sequence) :]
    return any(stem and post_guide.startswith(stem) for stem in stems)


def select_representative(
    candidates: Sequence[Mapping[str, object]],
    raw_records: Mapping[str, Mapping[str, str]],
    stems: Sequence[str],
) -> Mapping[str, object]:
    """Choose the dominant exact construct, breaking ties by truth ID."""

    eligible = [
        row
        for row in candidates
        if int(row.get("mapq", -1)) == 255
        and int(row.get("nm", -1)) == 0
        and as_bool(row.get("target_concordant"))
        and as_bool(row.get("production_filter_removed"))
        and int(row.get("tso_offset", -1)) == 0
    ]
    if not eligible:
        raise ValueError("No representative candidates passed alignment/truth gates")

    patterns = Counter((str(row["guide_query"]), str(row["cigar"])) for row in eligible)
    ordered_patterns = sorted(patterns, key=lambda item: (-patterns[item], item[0], item[1]))
    for guide, cigar in ordered_patterns:
        rows = sorted(
            (
                row
                for row in eligible
                if str(row["guide_query"]) == guide and str(row["cigar"]) == cigar
            ),
            key=lambda row: str(row["truth_id"]),
        )
        for row in rows:
            original_id = str(row["original_read_id"])
            raw = raw_records.get(original_id)
            if raw and is_exact_construct(raw["sequence"], guide, stems):
                return row
    raise ValueError("No dominant-pattern candidate has an exact TSO-guide-stem construct")


def read_fastq_records(path: Path) -> Dict[str, Dict[str, str]]:
    """Read a compact truth FASTQ, keyed by the truth ID before ``|``."""

    records: Dict[str, Dict[str, str]] = {}
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline().rstrip("\r\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n").upper()
            plus = handle.readline().rstrip("\r\n")
            quality = handle.readline().rstrip("\r\n")
            if not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}: {header}")
            read_id = header.split()[0].lstrip("@").split("|", 1)[0]
            if read_id in records:
                raise ValueError(f"Duplicate truth ID in {path}: {read_id}")
            records[read_id] = {
                "header": header,
                "sequence": sequence,
                "quality": quality,
            }
    return records


def extract_raw_fastq_records(
    path: Path,
    wanted_ids: Sequence[str],
    threads: int,
    pigz: str,
) -> Dict[str, Dict[str, str]]:
    """Stream a large gzip FASTQ once and retain the requested records."""

    wanted = set(wanted_ids)
    found: Dict[str, Dict[str, str]] = {}
    process = subprocess.Popen(
        [pigz, "-dc", "-p", str(max(1, threads)), str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    while True:
        header = process.stdout.readline()
        if not header:
            break
        sequence = process.stdout.readline()
        plus = process.stdout.readline()
        quality = process.stdout.readline()
        if not sequence or not plus or not quality:
            process.kill()
            raise ValueError(f"Truncated FASTQ while reading {path}")
        read_id = header.split(maxsplit=1)[0].lstrip("@").rstrip("\r\n")
        if read_id.endswith(("/1", "/2")):
            read_id = read_id[:-2]
        if read_id in wanted:
            found[read_id] = {
                "header": header.rstrip("\r\n"),
                "sequence": sequence.rstrip("\r\n").upper(),
                "quality": quality.rstrip("\r\n"),
            }
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"pigz failed for {path}: {stderr.strip()}")
    missing = sorted(wanted - set(found))
    if missing:
        preview = ", ".join(missing[:5])
        raise RuntimeError(
            f"Raw FASTQ scan found {len(found)}/{len(wanted)} requested IDs; "
            f"missing examples: {preview}"
        )
    return found


def write_fastq(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for row in rows:
            handle.write(
                f"{row['header']}\n{row['sequence']}\n+\n{row['quality']}\n"
            )


def write_deterministic_gzip_text(path: Path, text: str) -> None:
    """Write reproducible gzip bytes without a timestamp or source filename."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode("utf-8"), compresslevel=9, mtime=0))


def write_deterministic_dataframe(frame: pd.DataFrame, path: Path) -> None:
    text = frame.to_csv(sep="\t", index=False, lineterminator="\n")
    if path.suffix == ".gz":
        write_deterministic_gzip_text(path, text)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)


def write_deterministic_fastq(
    path: Path, rows: Sequence[Mapping[str, str]]
) -> None:
    text = "".join(
        f"{row['header']}\n{row['sequence']}\n+\n{row['quality']}\n"
        for row in rows
    )
    write_deterministic_gzip_text(path, text)


def parse_gtf_attributes(value: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for field in value.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field or " " not in field:
            continue
        key, raw = field.split(" ", 1)
        attrs[key] = raw.strip().strip('"')
    return attrs


def collect_annotation_rows(genome_dir: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    by_id = {target.gene_id: target for target in TARGETS}
    for feature, filename in (("gene", "genes.gtf.gz"), ("exon", "exons.gtf.gz")):
        path = genome_dir / "process" / filename
        if not path.exists():
            raise FileNotFoundError(path)
        unique = set()
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != feature:
                    continue
                attrs = parse_gtf_attributes(fields[8])
                gene_id = attrs.get("gene_id", "")
                target = by_id.get(gene_id)
                if target is None or fields[0] != target.star_contig:
                    continue
                start, end = int(fields[3]), int(fields[4])
                if end < target.start_1 or start > target.end_1:
                    continue
                key = (target.gene, feature, start, end, fields[6])
                if key in unique:
                    continue
                unique.add(key)
                rows.append(
                    {
                        "gene": target.gene,
                        "gene_id": gene_id,
                        "contig": fields[0],
                        "feature": feature,
                        "start_1": start,
                        "end_1": end,
                        "strand": fields[6],
                    }
                )
    return sorted(rows, key=lambda row: (row["gene"], row["feature"], row["start_1"], row["end_1"]))


def resolve_reference_fasta(genome_dir: Path) -> Path:
    manifest_path = genome_dir / "process" / "mkref_def.json"
    with manifest_path.open() as handle:
        payload = json.load(handle)
    raw = payload["mkref"]["input_files"][0]["fasta"]
    path = Path(str(raw).split(",", 1)[0])
    if not path.exists() or not Path(f"{path}.fai").exists():
        raise FileNotFoundError(f"Indexed reference FASTA is unavailable: {path}")
    return path


def run_checked(command: Sequence[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr.strip()}"
        )


def extract_locus_bam(
    source_bam: Path,
    bed_path: Path,
    output_bam: Path,
    samtools: str,
    threads: int,
) -> None:
    """Stream an unindexed BAM and retain primary MAPQ-255 locus overlaps."""

    unsorted = output_bam.with_name(f"{output_bam.stem}.unsorted.bam")
    first = subprocess.Popen(
        [
            samtools,
            "view",
            "-@",
            str(max(1, threads)),
            "-h",
            "-F",
            "2308",
            "-q",
            "255",
            str(source_bam),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert first.stdout is not None
    second = subprocess.run(
        [
            samtools,
            "view",
            "-@",
            str(max(1, threads)),
            "-b",
            "-L",
            str(bed_path),
            "-o",
            str(unsorted),
            "-",
        ],
        stdin=first.stdout,
        capture_output=True,
    )
    first.stdout.close()
    first_stderr = first.stderr.read().decode("utf-8", errors="replace") if first.stderr else ""
    first_code = first.wait()
    if first_code != 0 or second.returncode != 0:
        raise RuntimeError(
            "samtools streaming extraction failed:\n"
            f"source view: {first_stderr.strip()}\n"
            f"locus view: {second.stderr.decode('utf-8', errors='replace').strip()}"
        )
    run_checked(
        [
            samtools,
            "sort",
            "-@",
            str(max(1, threads)),
            "-o",
            str(output_bam),
            str(unsorted),
        ]
    )
    run_checked([samtools, "index", "-@", str(max(1, threads)), str(output_bam)])
    unsorted.unlink(missing_ok=True)


def write_dataframe(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def validate_equal(label: str, observed: object, expected: object) -> None:
    if observed != expected:
        raise RuntimeError(f"{label}: observed {observed!r}, expected {expected!r}")


def load_yaml(path: Path) -> MutableMapping[str, object]:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def extract_command(args: argparse.Namespace) -> None:
    try:
        import pysam
        import splitpipe  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "extract requires the production splitpipe runtime with pysam"
        ) from exc

    from qc_config_loader import load_qc_config
    from qc_gex_grna_contamination import GeneAnnotator, make_model, reverse_complement

    config_path = Path(args.config).resolve()
    cfg = load_qc_config(str(config_path))
    qc_root = Path(args.qc_root or cfg["gex_grna_contamination_qc"]["output_dir"])
    output_dir = Path(args.output_dir).resolve()
    scratch_dir = Path(args.scratch_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)

    library_dir = qc_root / "C2"
    manifest_path = library_dir / "alignment_input" / "ground_truth_manifest.tsv.gz"
    full_fastq = library_dir / "alignment_input" / "full" / "C2_ground_truth_R1_full.fastq.gz"
    short_fastq = library_dir / "alignment_input" / "short54" / "C2_ground_truth_R1_54.fastq.gz"
    full_bam = library_dir / "alignment_star" / "full" / "C2_Aligned.out.bam"
    short_dir = library_dir / "alignment_star" / "short54"
    short_bam = short_dir / "C2_Aligned.out.bam"
    short_summary_path = short_dir / "alignment_impact_summary.yaml"
    full_summary_path = library_dir / "alignment_star" / "full" / "alignment_impact_summary.yaml"
    benchmark_path = library_dir / "benchmark_summary.yaml"
    molecule_path = short_dir / "false_gene_assigned_molecules_by_gene.csv"
    for path in (
        manifest_path,
        full_fastq,
        short_fastq,
        full_bam,
        short_bam,
        short_summary_path,
        full_summary_path,
        benchmark_path,
        molecule_path,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    inferred_target_json = config_path.parent.parent / "target_annotation" / "short54_target_concordance.json"
    target_json_path = Path(args.target_annotation_json).resolve() if args.target_annotation_json else inferred_target_json
    if not target_json_path.exists():
        raise FileNotFoundError(target_json_path)

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    validate_equal("truth panel input", len(manifest), EXPECTED_GLOBAL["input_reads"])
    manifest_by_truth = manifest.set_index("truth_id").to_dict("index")

    genome_dir = Path(cfg["paths"]["genome_dir"])
    gene_info_path = genome_dir / "gene_info.json"
    annotator = GeneAnnotator(str(gene_info_path))
    model = make_model(cfg, str(config_path))

    controlled_rows: List[Dict[str, object]] = []
    mapped_primary = unique = false_assigned = target_assigned = 0
    with pysam.AlignmentFile(str(short_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                continue
            mapped_primary += 1
            unique += int(read.mapping_quality == 255)
            truth = parse_truth_qname(read.query_name)
            truth_id = truth["truth_id"]
            meta = manifest_by_truth.get(truth_id)
            if meta is None:
                raise RuntimeError(f"Controlled BAM truth ID missing from manifest: {truth_id}")
            if read.mapping_quality == 255:
                gene_id, gene_name, exonic_fraction = annotator.assign(read)
                gene_name = normalize_gene_name(gene_name)
            else:
                gene_id, gene_name, exonic_fraction = "", "", 0.0
            false_assigned += int(bool(gene_id))
            expected_targets = {value for value in str(meta["target_genes"]).split(",") if value}
            target_concordant = bool(gene_name and gene_name in expected_targets)
            target_assigned += int(target_concordant)
            blocks = [(start + 1, end) for start, end in read.get_blocks()]
            row = {
                "truth_id": truth_id,
                "original_read_id": meta["original_read_id"],
                "cell": meta["cell"],
                "barcode_indices": meta["barcode_indices"],
                "umi": meta["umi"],
                "truth_tier": meta["truth_tier"],
                "guide_query": meta["guide_query"],
                "guide_names": meta["guide_names"],
                "target_genes": meta["target_genes"],
                "tso_offset": int(meta["tso_offset"]),
                "production_filter_removed": str(meta["production_filter_removed"]).lower() == "true",
                "contig": read.reference_name,
                "start_1": read.reference_start + 1,
                "end_1": read.reference_end,
                "blocks": serialize_blocks(blocks),
                "cigar": read.cigarstring,
                "strand": "-" if read.is_reverse else "+",
                "mapq": int(read.mapping_quality),
                "nh": int(read.get_tag("NH")) if read.has_tag("NH") else -1,
                "nm": int(read.get_tag("nM")) if read.has_tag("nM") else (
                    int(read.get_tag("NM")) if read.has_tag("NM") else -1
                ),
                "star_sequence": read.query_sequence or "",
                "assigned_gene_id": gene_id,
                "assigned_gene_name": gene_name,
                "exonic_fraction": float(exonic_fraction),
                "annotation_region": "exonic" if exonic_fraction >= 0.5 else ("intronic" if gene_id else "unassigned"),
                "target_concordant": target_concordant,
            }
            controlled_rows.append(row)

    validate_equal("short54 mapped primary", mapped_primary, EXPECTED_GLOBAL["short_mapped_primary_reads"])
    validate_equal("short54 unique", unique, EXPECTED_GLOBAL["short_unique_reads"])
    validate_equal("short54 false assigned", false_assigned, EXPECTED_GLOBAL["short_false_gene_assigned_reads"])
    validate_equal("short54 target assigned", target_assigned, EXPECTED_GLOBAL["short_guide_target_assigned_reads"])
    validate_equal(
        "short54 multimapping primary",
        mapped_primary - unique,
        EXPECTED_GLOBAL["short_multimapping_primary_reads"],
    )

    full_mapped = 0
    with pysam.AlignmentFile(str(full_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
                full_mapped += 1
    validate_equal("full mapped primary", full_mapped, EXPECTED_GLOBAL["full_mapped_primary_reads"])

    target_payload = json.loads(target_json_path.read_text())
    if not target_payload.get("complete"):
        raise RuntimeError(f"Corrected target annotation is incomplete: {target_json_path}")
    for key, expected_key in (
        ("mapped_primary_reads", "short_mapped_primary_reads"),
        ("unique_reads", "short_unique_reads"),
        ("gene_assigned_reads", "short_false_gene_assigned_reads"),
        ("guide_target_assigned_reads", "short_guide_target_assigned_reads"),
    ):
        validate_equal(f"corrected target annotation {key}", int(target_payload["counts"][key]), EXPECTED_GLOBAL[expected_key])

    controlled = pd.DataFrame(controlled_rows)
    for target in TARGETS:
        count = int(
            (
                (controlled["assigned_gene_name"] == target.gene)
                & controlled["target_concordant"]
            ).sum()
        )
        validate_equal(f"{target.gene} target-concordant reads", count, target.expected_target_reads)

    molecule_frame = pd.read_csv(molecule_path)
    molecule_frame["gene_name"] = molecule_frame["gene_name"].map(normalize_gene_name)
    removed_molecules = {}
    for target in TARGETS:
        values = molecule_frame.loc[
            (molecule_frame["gene_name"] == target.gene)
            & (molecule_frame["production_filter_stratum"] == "removed"),
            "false_gene_assigned_molecules",
        ]
        observed = int(values.sum())
        validate_equal(f"{target.gene} removed false molecules", observed, target.expected_removed_false_molecules)
        removed_molecules[target.gene] = observed

    candidate_rows = []
    for target in TARGETS:
        subset = controlled.loc[
            (controlled["assigned_gene_name"] == target.gene)
            & controlled["target_concordant"]
            & (controlled["mapq"] == 255)
            & (controlled["nm"] == 0)
            & controlled["production_filter_removed"]
            & (controlled["tso_offset"] == 0)
        ]
        candidate_rows.extend(subset.to_dict("records"))
    wanted_original_ids = sorted({str(row["original_read_id"]) for row in candidate_rows})

    benchmark = load_yaml(benchmark_path)
    raw_r1 = Path(benchmark["input_fastq_pairs"][0]["r1"])
    pigz = args.pigz or shutil.which("pigz")
    if not pigz:
        raise RuntimeError("pigz is required for the one-pass raw FASTQ extraction")
    raw_records = extract_raw_fastq_records(raw_r1, wanted_original_ids, args.threads, pigz)

    stems = [
        str(cfg["grna_structure"]["fixed_seq_2"]).upper(),
        str(cfg["grna_structure"]["fixed_seq_3"]).upper(),
    ]
    full_records = read_fastq_records(full_fastq)
    short_records = read_fastq_records(short_fastq)
    reference_fasta = Path(args.reference_fasta).resolve() if args.reference_fasta else resolve_reference_fasta(genome_dir)
    fasta = pysam.FastaFile(str(reference_fasta))

    representative_rows: List[Dict[str, object]] = []
    representative_fastq_rows: List[Dict[str, str]] = []
    for target in TARGETS:
        target_candidates = [row for row in candidate_rows if row["assigned_gene_name"] == target.gene]
        selected = dict(select_representative(target_candidates, raw_records, stems))
        truth_id = str(selected["truth_id"])
        original_id = str(selected["original_read_id"])
        raw = raw_records[original_id]
        short54 = raw["sequence"][:54]
        processed_short, trim_index = splitpipe_trim(short54)
        processed_full, full_trim_index = splitpipe_trim(raw["sequence"])
        validate_equal(f"{target.gene} short54 FASTQ", processed_short, short_records[truth_id]["sequence"])
        validate_equal(f"{target.gene} full FASTQ", processed_full, full_records[truth_id]["sequence"])
        validate_equal(f"{target.gene} short processed length", len(processed_short), 24)
        validate_equal(f"{target.gene} full processed length", len(processed_full), 120)
        validate_equal(f"{target.gene} TSO anchor index", trim_index, 13)
        validate_equal(f"{target.gene} full TSO anchor index", full_trim_index, 13)

        aligned_pairs = []
        query_pos = 0
        ref_pos = int(selected["start_1"])
        reference_aligned = []
        query_aligned = []
        match_symbols = []
        for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", str(selected["cigar"])):
            length = int(length_text)
            if op in "M=X":
                ref_seq = fasta.fetch(target.fasta_contig, ref_pos - 1, ref_pos - 1 + length).upper()
                query_seq = str(selected["star_sequence"])[query_pos : query_pos + length]
                reference_aligned.append(ref_seq)
                query_aligned.append(query_seq)
                match_symbols.append("".join("|" if left == right else "x" for left, right in zip(query_seq, ref_seq)))
                aligned_pairs.extend(zip(query_seq, ref_seq))
                query_pos += length
                ref_pos += length
            elif op == "I":
                query_pos += length
            elif op in "DN":
                ref_pos += length
            elif op == "S":
                query_pos += length
        mismatches = sum(left != right for left, right in aligned_pairs)
        validate_equal(f"{target.gene} reference mismatches", mismatches, 0)
        validate_equal(f"{target.gene} NM", int(selected["nm"]), 0)

        selected.update(
            {
                "raw_r1": raw["sequence"],
                "synthetic54_r1": short54,
                "processed_short54_r1": processed_short,
                "processed_full_r1": processed_full,
                "raw_tso": raw["sequence"][: len(TSO)],
                "raw_guide": raw["sequence"][len(TSO) : len(TSO) + len(str(selected["guide_query"]))],
                "raw_post_guide": raw["sequence"][len(TSO) + len(str(selected["guide_query"])) : 54],
                "reference_aligned": "".join(reference_aligned),
                "query_aligned": "".join(query_aligned),
                "match_symbols": "".join(match_symbols),
                "full_arm_mapped": False,
            }
        )
        representative_rows.append(selected)
        representative_fastq_rows.append(
            {
                "header": f"@{truth_id}|original={original_id}|target={target.gene}",
                "sequence": raw["sequence"],
                "quality": raw["quality"],
            }
        )
    fasta.close()

    bed_path = scratch_dir / "c2_case_loci.bed"
    with bed_path.open("w") as handle:
        for target in TARGETS:
            handle.write(f"{target.star_contig}\t{target.start_1 - 1}\t{target.end_1}\t{target.gene}\n")
    source_gex_bam = Path(args.gex_bam).resolve() if args.gex_bam else Path(cfg["paths"]["output_dir"]) / "gex_C2" / "process" / "barcode_headAligned_anno.bam"
    if not source_gex_bam.exists():
        raise FileNotFoundError(source_gex_bam)
    samtools = args.samtools or shutil.which("samtools")
    if not samtools:
        raise RuntimeError("samtools is required for production-BAM extraction")
    mini_bam = scratch_dir / "c2_case_loci.sorted.bam"
    extract_locus_bam(source_gex_bam, bed_path, mini_bam, samtools, args.threads)

    target_by_contig = {target.star_contig: target for target in TARGETS}
    gex_rows: List[Dict[str, object]] = []
    excluded_structural = Counter()
    excluded_span_only = Counter()
    locus_primary_mapq255 = Counter()
    with pysam.AlignmentFile(str(mini_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                raise RuntimeError("Mini-BAM contains a non-primary mapped record")
            if read.mapping_quality != 255:
                raise RuntimeError("Mini-BAM contains a non-MAPQ-255 record")
            target = target_by_contig.get(read.reference_name)
            if target is None:
                continue
            blocks = [(start + 1, end) for start, end in read.get_blocks()]
            if not intervals_overlap(blocks, target.start_1, target.end_1):
                # ``samtools view -L`` uses the outer reference span, so a
                # splice can cross the BED window without an aligned block in
                # it. Such span-only candidates are not locus reads here.
                excluded_span_only[target.gene] += 1
                continue
            locus_primary_mapq255[target.gene] += 1
            gx = read.get_tag("GX") if read.has_tag("GX") else ""
            if gx != target.gene_id:
                continue
            sequence = read.query_sequence or ""
            classifier_sequence = reverse_complement(sequence) if read.is_reverse else sequence
            evidence = model.classify_post_splitpipe(classifier_sequence)
            if evidence.structure_truth:
                excluded_structural[target.gene] += 1
                continue
            gene_name = normalize_gene_name(read.get_tag("GN") if read.has_tag("GN") else "")
            if gene_name != target.gene:
                raise RuntimeError(
                    f"GX/GN mismatch for {read.query_name}: {gx}/{gene_name}, expected {target.gene_id}/{target.gene}"
                )
            gex_rows.append(
                {
                    "gene": target.gene,
                    "gene_id": gx,
                    "read_id": read.query_name,
                    "contig": read.reference_name,
                    "start_1": read.reference_start + 1,
                    "end_1": read.reference_end,
                    "blocks": serialize_blocks(blocks),
                    "cigar": read.cigarstring,
                    "strand": "-" if read.is_reverse else "+",
                    "mapq": int(read.mapping_quality),
                    "nm": int(read.get_tag("nM")) if read.has_tag("nM") else (
                        int(read.get_tag("NM")) if read.has_tag("NM") else -1
                    ),
                    "nh": int(read.get_tag("NH")) if read.has_tag("NH") else -1,
                    "gx": gx,
                    "gn": gene_name,
                    "cell": read.get_tag("CB") if read.has_tag("CB") else "",
                    "barcode_indices": read.get_tag("pB") if read.has_tag("pB") else "",
                    "umi": read.get_tag("pN") if read.has_tag("pN") else "",
                    "query_length": int(read.query_length or 0),
                    "query_sequence": sequence,
                    "post_splitpipe_structure_truth": False,
                }
            )

    annotation_rows = collect_annotation_rows(genome_dir)
    fasta = pysam.FastaFile(str(reference_fasta))
    reference_rows = []
    for target in TARGETS:
        sequence = fasta.fetch(target.fasta_contig, target.start_1 - 1, target.end_1).upper()
        validate_equal(f"{target.gene} reference window length", len(sequence), target.end_1 - target.start_1 + 1)
        reference_rows.append(
            {
                "gene": target.gene,
                "star_contig": target.star_contig,
                "fasta_contig": target.fasta_contig,
                "start_1": target.start_1,
                "end_1": target.end_1,
                "sequence": sequence,
            }
        )
    fasta.close()

    gex = pd.DataFrame(gex_rows)
    locus_rows = []
    for target in TARGETS:
        target_controlled = controlled.loc[
            (controlled["assigned_gene_name"] == target.gene)
            & controlled["target_concordant"]
        ].copy()
        pattern_counts = (
            target_controlled.groupby(["guide_query", "start_1", "cigar"], dropna=False)
            .size()
            .reset_index(name="read_count")
            .sort_values(["read_count", "guide_query", "start_1", "cigar"], ascending=[False, True, True, True])
        )
        pattern_text = "; ".join(
            f"{row.guide_query}@{target.star_contig}:{int(row.start_1)} {row.cigar} (n={int(row.read_count)})"
            for row in pattern_counts.itertuples(index=False)
        )
        regions = sorted(set(target_controlled["annotation_region"]))
        validate_equal(f"{target.gene} annotation class", regions, [target.expected_region])
        locus_rows.append(
            {
                "gene": target.gene,
                "gene_id": target.gene_id,
                "contig": target.star_contig,
                "window_start_1": target.start_1,
                "window_end_1": target.end_1,
                "annotation_region": target.expected_region,
                "target_concordant_reads": len(target_controlled),
                "removed_false_molecules": removed_molecules[target.gene],
                "guide_alignment_patterns": pattern_text,
                "putative_genuine_gex_reads": int((gex["gene"] == target.gene).sum()) if not gex.empty else 0,
                "all_primary_mapq255_locus_reads": int(locus_primary_mapq255[target.gene]),
                "excluded_span_only_locus_candidates": int(excluded_span_only[target.gene]),
                "excluded_structural_gex_bam_reads": int(excluded_structural[target.gene]),
            }
        )

    benchmark_raw = benchmark["raw_read_prevalence"]
    lower = benchmark["prealignment_molecule_prevalence"]["fractions"]["cell_pb_umi_insert"]["contaminant_of_all"]["fraction"] * 1_000_000
    upper = benchmark["prealignment_molecule_prevalence"]["fractions"]["cell_pb_umi"]["contaminant_of_all"]["fraction"] * 1_000_000
    global_summary = {
        "library": "C2",
        "raw_reads": int(benchmark["records_evaluated"]),
        "raw_structural_contaminants": int(benchmark_raw["structure_contaminant_reads"]),
        "raw_structural_contaminant_cpm": float(benchmark_raw["contaminant_fraction_of_all_reads"] * 1_000_000),
        "prealignment_family_cpm_low": float(lower),
        "prealignment_family_cpm_high": float(upper),
        **EXPECTED_GLOBAL,
        "guide_target_fraction_of_input": EXPECTED_GLOBAL["short_guide_target_assigned_reads"] / EXPECTED_GLOBAL["input_reads"],
        "truth_panel_token": benchmark["alignment_ground_truth_input"]["truth_panel_token"],
    }

    write_dataframe(controlled.sort_values(["contig", "start_1", "truth_id"]), output_dir / "controlled_short54_alignments.tsv.gz")
    write_dataframe(pd.DataFrame(representative_rows), output_dir / "representative_reads.tsv")
    write_fastq(output_dir / "representative_raw_R1.fastq.gz", representative_fastq_rows)
    write_dataframe(
        gex.sort_values(["gene", "start_1", "end_1", "read_id"]) if not gex.empty else pd.DataFrame(columns=[
            "gene", "gene_id", "read_id", "contig", "start_1", "end_1", "blocks", "cigar", "strand", "mapq", "nm", "nh", "gx", "gn", "cell", "barcode_indices", "umi", "query_length", "query_sequence", "post_splitpipe_structure_truth"
        ]),
        output_dir / "putative_genuine_gex_reads.tsv.gz",
    )
    write_dataframe(pd.DataFrame(annotation_rows), output_dir / "locus_annotation.tsv")
    write_dataframe(pd.DataFrame(reference_rows), output_dir / "reference_windows.tsv.gz")
    write_dataframe(pd.DataFrame(locus_rows), output_dir / "locus_summary.tsv")
    (output_dir / "global_summary.json").write_text(json.dumps(global_summary, indent=2, sort_keys=True) + "\n")

    provenance = {
        "complete": True,
        "job_id": os.environ.get("JOB_ID", "interactive"),
        "threads": int(args.threads),
        "config": str(config_path),
        "qc_root": str(qc_root),
        "truth_manifest": str(manifest_path),
        "controlled_short54_bam": str(short_bam),
        "controlled_full_bam": str(full_bam),
        "corrected_target_annotation": str(target_json_path),
        "raw_r1": str(raw_r1),
        "source_gex_bam": str(source_gex_bam),
        "parse_genome_dir": str(genome_dir),
        "reference_fasta": str(reference_fasta),
        "truth_panel_token": global_summary["truth_panel_token"],
        "output_files": {},
    }
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "extraction.complete.json":
            provenance["output_files"][path.name] = {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    (output_dir / "extraction.complete.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    mini_bam.unlink(missing_ok=True)
    Path(f"{mini_bam}.bai").unlink(missing_ok=True)
    print(f"Case-study extraction complete: {output_dir}")


def load_bundle(data_dir: Path):
    required = [
        "extraction.complete.json",
        "global_summary.json",
        "controlled_short54_alignments.tsv.gz",
        "representative_reads.tsv",
        "representative_raw_R1.fastq.gz",
        "putative_genuine_gex_reads.tsv.gz",
        "locus_annotation.tsv",
        "reference_windows.tsv.gz",
        "locus_summary.tsv",
    ]
    for name in required:
        if not (data_dir / name).exists():
            raise FileNotFoundError(data_dir / name)
    provenance = json.loads((data_dir / "extraction.complete.json").read_text())
    if not provenance.get("complete"):
        raise RuntimeError("Extraction completion manifest is not complete")
    output_files = provenance.get("output_files")
    if not isinstance(output_files, dict):
        raise RuntimeError("Extraction manifest lacks hashed output metadata")
    unhashed = sorted(set(required[1:]) - set(output_files))
    if unhashed:
        raise RuntimeError(
            "Extraction manifest lacks hashes for required files: "
            + ", ".join(unhashed)
        )
    for name, metadata in output_files.items():
        path = data_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        validate_equal(f"SHA256 {name}", sha256_file(path), metadata["sha256"])
    return {
        "provenance": provenance,
        "global": json.loads((data_dir / "global_summary.json").read_text()),
        "controlled": pd.read_csv(data_dir / "controlled_short54_alignments.tsv.gz", sep="\t"),
        "representatives": pd.read_csv(data_dir / "representative_reads.tsv", sep="\t"),
        "gex": pd.read_csv(data_dir / "putative_genuine_gex_reads.tsv.gz", sep="\t"),
        "annotation": pd.read_csv(data_dir / "locus_annotation.tsv", sep="\t"),
        "reference": pd.read_csv(data_dir / "reference_windows.tsv.gz", sep="\t"),
        "loci": pd.read_csv(data_dir / "locus_summary.tsv", sep="\t"),
    }


def _published_columns(
    frame: pd.DataFrame, columns: Sequence[str], label: str
) -> pd.DataFrame:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} lacks publishable columns: {', '.join(missing)}")
    return frame.loc[:, list(columns)].copy()


def publish_sanitized_bundle(source_dir: Path, output_dir: Path) -> Path:
    """Publish the scientific minimum needed to reproduce the case study.

    Direct instrument identifiers, cell/barcode/UMI fields, unused native-GEX
    sequences, and source paths are deliberately omitted. Stable sequential
    evidence IDs preserve deterministic sorting without retaining a mapping to
    the source identifiers.
    """

    source_dir = source_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"published evidence directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    bundle = load_bundle(source_dir)
    controlled = bundle["controlled"].copy()
    representatives = bundle["representatives"].copy()
    gex = bundle["gex"].copy()

    controlled_ids = {
        value: f"C2CTL{index:06d}"
        for index, value in enumerate(
            sorted(controlled["truth_id"].astype(str).unique()), start=1
        )
    }
    controlled["truth_id"] = controlled["truth_id"].astype(str).map(controlled_ids)
    controlled = _published_columns(
        controlled, PUBLISHED_CONTROLLED_COLUMNS, "controlled alignments"
    ).sort_values(["contig", "start_1", "truth_id"], kind="stable")

    representative_source_ids = representatives["truth_id"].astype(str).tolist()
    representative_ids = {
        value: f"C2REP{index:03d}"
        for index, value in enumerate(sorted(representative_source_ids), start=1)
    }
    representatives["truth_id"] = (
        representatives["truth_id"].astype(str).map(representative_ids)
    )
    representatives = _published_columns(
        representatives,
        PUBLISHED_REPRESENTATIVE_COLUMNS,
        "representative reads",
    ).sort_values(["assigned_gene_name", "truth_id"], kind="stable")

    gex = gex.sort_values(
        ["gene", "start_1", "end_1", "read_id"], kind="stable"
    ).reset_index(drop=True)
    gex["read_id"] = [f"C2GEX{index:06d}" for index in range(1, len(gex) + 1)]
    gex = _published_columns(gex, PUBLISHED_GEX_COLUMNS, "putative GEX reads")

    raw_fastq = read_fastq_records(source_dir / "representative_raw_R1.fastq.gz")
    published_fastq = []
    for source_id in sorted(representative_source_ids):
        if source_id not in raw_fastq:
            raise ValueError(f"representative FASTQ lacks truth ID {source_id}")
        row = representatives.loc[
            representatives["truth_id"] == representative_ids[source_id]
        ].iloc[0]
        published_fastq.append(
            {
                "header": (
                    f"@{representative_ids[source_id]}|"
                    f"target={row['assigned_gene_name']}"
                ),
                "sequence": raw_fastq[source_id]["sequence"],
                "quality": raw_fastq[source_id]["quality"],
            }
        )

    frames = {
        "controlled_short54_alignments.tsv.gz": controlled,
        "representative_reads.tsv": representatives,
        "putative_genuine_gex_reads.tsv.gz": gex,
        "locus_annotation.tsv": bundle["annotation"].sort_values(
            ["gene", "feature", "start_1", "end_1"], kind="stable"
        ),
        "reference_windows.tsv.gz": bundle["reference"].sort_values(
            ["gene", "start_1", "end_1"], kind="stable"
        ),
        "locus_summary.tsv": bundle["loci"].sort_values("gene", kind="stable"),
    }
    for name, frame in frames.items():
        write_deterministic_dataframe(frame, output_dir / name)
    write_deterministic_fastq(
        output_dir / "representative_raw_R1.fastq.gz", published_fastq
    )
    (output_dir / "global_summary.json").write_text(
        json.dumps(bundle["global"], indent=2, sort_keys=True) + "\n"
    )

    source_manifest = source_dir / "extraction.complete.json"
    provenance = bundle["provenance"]
    manifest = {
        "schema_version": PUBLISHED_BUNDLE_SCHEMA,
        "complete": True,
        "sanitized": True,
        "library": "C2",
        "job_id": provenance.get("job_id", "unknown"),
        "truth_panel_token": bundle["global"]["truth_panel_token"],
        "source_extraction_manifest_sha256": sha256_file(source_manifest),
        "privacy": {
            "direct_read_identifiers": "pseudonymized_without_mapping",
            "cell_barcode_umi_fields": "removed",
            "unused_native_gex_query_sequences": "removed",
            "absolute_source_paths": "removed",
            "representative_scientific_sequences": "retained",
        },
        "output_files": {},
    }
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "extraction.complete.json":
            manifest["output_files"][path.name] = {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    manifest_path = output_dir / "extraction.complete.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest_path


def publish_command(args: argparse.Namespace) -> None:
    manifest = publish_sanitized_bundle(
        Path(args.data_dir), Path(args.output_dir)
    )
    print(f"Sanitized case-study evidence published: {manifest}")


def setup_matplotlib():
    cache_dir = Path(tempfile.gettempdir()) / "csm_c2_case_study_matplotlib"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.hashsalt"] = "c2-gex-grna-misalignment"
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    return plt, Rectangle


def save_figure(fig, asset_dir: Path, stem: str) -> None:
    svg_path = asset_dir / f"{stem}.svg"
    fig.savefig(
        svg_path,
        bbox_inches="tight",
        metadata={"Creator": "c2_gex_grna_case_study.py", "Date": None},
    )
    svg_text = svg_path.read_text()
    svg_path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n")
    fig.savefig(
        asset_dir / f"{stem}.png",
        dpi=180,
        bbox_inches="tight",
        metadata={"Software": "c2_gex_grna_case_study.py"},
    )


def aligned_base_components(length: int, guide_length: int) -> List[str]:
    """Annotate aligned positions by their source guide/stem component."""

    return ["guide" if index < guide_length else "stem" for index in range(length)]


def add_colored_alignment(
    axis, query: str, symbols: str, reference: str, guide_length: int
) -> None:
    """Color query and reference letters by the query's source component."""

    from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker

    label_props = {"fontsize": 9, "fontfamily": "monospace", "color": "#222222"}

    components = aligned_base_components(len(query), guide_length)

    def sequence_row(label: str, sequence: str):
        cells = [TextArea(f"{label:<8}", textprops=label_props)]
        cells.extend(
            TextArea(
                base,
                textprops={
                    "fontsize": 9,
                    "fontfamily": "monospace",
                    "fontweight": "bold",
                    "color": COMPONENT_COLORS[component],
                },
            )
            for base, component in zip(sequence, components)
        )
        return HPacker(children=cells, align="center", pad=0, sep=0)

    match_cells = [TextArea(" " * 8, textprops=label_props)]
    match_cells.extend(
        TextArea(
            symbol,
            textprops={
                "fontsize": 9,
                "fontfamily": "monospace",
                "color": "#2E8B57" if symbol == "|" else "#D62728",
            },
        )
        for symbol in symbols
    )
    rows = [
        sequence_row("Query", query),
        HPacker(children=match_cells, align="center", pad=0, sep=0),
        sequence_row("GRCh38", reference),
    ]
    block = VPacker(children=rows, align="left", pad=0, sep=-1)
    artist = AnnotationBbox(
        block,
        (0.0, 0.055),
        xycoords=axis.transAxes,
        box_alignment=(0, 0),
        frameon=False,
        pad=0,
    )
    axis.add_artist(artist)


def render_transformation(bundle: Mapping[str, object], asset_dir: Path) -> None:
    plt, Rectangle = setup_matplotlib()
    representatives = bundle["representatives"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 8.0))
    colors = COMPONENT_COLORS
    for axis, row in zip(axes, representatives.sort_values("assigned_gene_name").itertuples(index=False)):
        guide_length = len(str(row.guide_query))
        stem_visible = 54 - len(TSO) - guide_length
        stages = [("Raw C2 R1", 150), ("Synthetic R1=54", 54), ("STAR input", len(str(row.processed_short54_r1)))]
        axis.set_xlim(0, 150)
        axis.set_ylim(-0.5, 3.6)
        y_positions = [2.8, 1.8, 0.8]
        for (label, length), y in zip(stages, y_positions):
            axis.text(-2, y + 0.18, label, ha="right", va="center", fontsize=9)
            if label == "Raw C2 R1":
                segments = [
                    (start, end - start, component)
                    for start, end, component in raw_component_segments(
                        str(row.raw_r1), str(row.guide_query)
                    )
                ]
            elif label == "Synthetic R1=54":
                segments = [
                    (0, len(TSO), "TSO"),
                    (len(TSO), guide_length, "guide"),
                    (len(TSO) + guide_length, stem_visible, "stem"),
                ]
            else:
                segments = [
                    (0, min(guide_length, length), "guide"),
                    (guide_length, max(0, length - guide_length), "stem"),
                ]
            for start, width, segment in segments:
                if width <= 0:
                    continue
                axis.add_patch(Rectangle((start, y), width, 0.38, facecolor=colors[segment], edgecolor="white", linewidth=0.8))
                if label == "Raw C2 R1" and width >= 7:
                    short_labels = {
                        "TSO": "TSO",
                        "guide": "guide",
                        "stem": "stem",
                        "ibar": "iBAR",
                        "tracr": "tracr",
                        "bridge": "bridge",
                        "r2_barcode_rc": "BC RC",
                        "r2_linker_rc": "link RC",
                    }
                    axis.text(
                        start + width / 2,
                        y + 0.19,
                        short_labels[segment],
                        ha="center",
                        va="center",
                        fontsize=5.5,
                        color="white" if segment not in {"stem", "bridge", "r2_linker_rc"} else "#222222",
                    )
            axis.text(length + 1, y + 0.18, f"{length} nt", va="center", fontsize=8)

        add_colored_alignment(
            axis,
            str(row.query_aligned),
            str(row.match_symbols),
            str(row.reference_aligned),
            guide_length,
        )
        axis.set_title(
            f"{row.assigned_gene_name}: {row.truth_id}\n"
            f"{row.contig}:{int(row.start_1):,} {row.cigar}, strand {row.strand}, "
            f"MAPQ {int(row.mapq)}, NM {int(row.nm)}",
            fontsize=11,
        )
        axis.set_xlabel("R1 base position before/after preprocessing")
        axis.set_yticks([])
        axis.spines[["left", "right", "top"]].set_visible(False)
        axis.set_xticks([0, 20, 30, 54, 100, 150])
    component_order = (
        "TSO",
        "guide",
        "stem",
        "ibar",
        "tracr",
        "bridge",
        "r2_barcode_rc",
        "r2_linker_rc",
    )
    handles = [Rectangle((0, 0), 1, 1, facecolor=colors[key]) for key in component_order]
    fig.legend(
        handles,
        [
            "Split_TSO",
            "guide",
            "CSM stem",
            "iBAR",
            "partial tracr",
            "post-tracr bridge",
            "R2 barcode RC",
            "R2 linker RC",
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=4,
        fontsize=8,
    )
    fig.suptitle(
        "How a full TSO–gRNA contaminant becomes a short STAR alignment",
        fontsize=14,
        y=0.98,
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.20, top=0.83, wspace=0.30)
    save_figure(fig, asset_dir, "read_transformation")
    plt.close(fig)


def coverage_vector(rows: pd.DataFrame, start_1: int, end_1: int) -> List[int]:
    coverage = [0] * (end_1 - start_1 + 1)
    for value in rows.get("blocks", pd.Series(dtype=str)).fillna(""):
        for block_start, block_end in deserialize_blocks(value):
            left = max(start_1, block_start)
            right = min(end_1, block_end)
            for position in range(left, right + 1):
                coverage[position - start_1] += 1
    return coverage


def render_locus(bundle: Mapping[str, object], target: TargetLocus, asset_dir: Path) -> None:
    plt, Rectangle = setup_matplotlib()
    controlled = bundle["controlled"]
    gex = bundle["gex"]
    annotation = bundle["annotation"]
    representatives = bundle["representatives"]
    locus_summary = bundle["loci"].set_index("gene").loc[target.gene]
    rep = representatives.loc[representatives["assigned_gene_name"] == target.gene].iloc[0]
    target_controlled = controlled.loc[
        (controlled["assigned_gene_name"] == target.gene)
        & controlled["target_concordant"].map(as_bool)
    ]
    target_gex = gex.loc[gex["gene"] == target.gene] if not gex.empty else gex
    target_annotation = annotation.loc[annotation["gene"] == target.gene]

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(14, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.6, 1.6, 2.4]},
        constrained_layout=True,
    )
    positions = list(range(target.start_1, target.end_1 + 1))

    gene_rows = target_annotation.loc[target_annotation["feature"] == "gene"]
    gene_start = int(gene_rows["start_1"].min()) if not gene_rows.empty else target.start_1
    gene_end = int(gene_rows["end_1"].max()) if not gene_rows.empty else target.end_1
    axes[0].plot([gene_start, gene_end], [0.5, 0.5], color="#444444", linewidth=1.5)
    exon_intervals = merge_intervals(
        zip(
            target_annotation.loc[target_annotation["feature"] == "exon", "start_1"],
            target_annotation.loc[target_annotation["feature"] == "exon", "end_1"],
        )
    )
    for start, end in exon_intervals:
        axes[0].add_patch(Rectangle((start, 0.25), end - start + 1, 0.5, facecolor="#555555"))
    axes[0].text(
        0.01,
        0.92,
        f"Parse GTF: {target.gene} ({target.expected_region} guide assignment)",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
    )
    axes[0].set_ylim(0, 1)
    axes[0].set_yticks([])

    contaminant_coverage = coverage_vector(target_controlled, target.start_1, target.end_1)
    axes[1].fill_between(positions, contaminant_coverage, step="mid", color="#D62728", alpha=0.8)
    axes[1].set_ylabel("Contaminant\ncoverage")
    axes[1].text(
        0.99,
        0.9,
        f"n={int(locus_summary.target_concordant_reads):,} target-concordant reads",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )

    gex_coverage = coverage_vector(target_gex, target.start_1, target.end_1)
    axes[2].fill_between(positions, gex_coverage, step="mid", color="#4C78A8", alpha=0.8)
    axes[2].set_ylabel("Putative GEX\ncoverage")
    axes[2].text(
        0.99,
        0.9,
        f"n={int(locus_summary.putative_genuine_gex_reads):,} primary MAPQ-255 target-GX reads",
        transform=axes[2].transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    if target_gex.empty:
        axes[2].text(0.5, 0.5, "No qualifying reads observed", transform=axes[2].transAxes, ha="center", va="center")

    display_rows = target_gex.sort_values(["start_1", "end_1", "read_id"]).head(30) if not target_gex.empty else target_gex
    for y, row in enumerate(display_rows.itertuples(index=False)):
        for start, end in deserialize_blocks(row.blocks):
            axes[3].add_patch(Rectangle((start, y - 0.32), end - start + 1, 0.64, facecolor="#4C78A8", edgecolor="none"))
    axes[3].set_ylim(-1, max(2, len(display_rows)))
    axes[3].invert_yaxis()
    axes[3].set_ylabel("Representative\nGEX reads")
    if len(target_gex) > len(display_rows):
        axes[3].text(0.99, 0.92, f"first {len(display_rows)} of {len(target_gex):,} shown", transform=axes[3].transAxes, ha="right", va="top", fontsize=8)
    if display_rows.empty:
        axes[3].text(0.5, 0.5, "Empty track", transform=axes[3].transAxes, ha="center", va="center")

    for axis in axes:
        axis.axvspan(int(rep.start_1), int(rep.end_1), color="#D62728", alpha=0.10)
        axis.spines[["right", "top"]].set_visible(False)
        axis.ticklabel_format(style="plain", axis="x", useOffset=False)
    axes[3].set_xlim(target.start_1, target.end_1)
    axes[3].set_xlabel(f"{target.star_contig} coordinate (GRCh38, 1-based)")
    fig.suptitle(
        f"{target.gene}: contaminant and putative genuine C2 GEX alignments on the same annotation\n"
        f"Representative {rep.truth_id}: {rep.query_aligned} → {rep.reference_aligned}; "
        f"{rep.cigar}, strand {rep.strand}, MAPQ {int(rep.mapq)}, NM {int(rep.nm)}",
        fontsize=13,
    )
    save_figure(fig, asset_dir, f"{target.gene.lower()}_locus")
    plt.close(fig)


def format_percent(value: float, digits: int = 2) -> str:
    return f"{100 * float(value):.{digits}f}%"


def summarize_positions(controlled: pd.DataFrame, target: TargetLocus) -> Tuple[str, str]:
    """Format per-position counts and dominant guide sequences for Markdown."""

    selected = controlled.loc[
        (controlled["assigned_gene_name"] == target.gene)
        & controlled["target_concordant"].map(as_bool)
    ]
    position_labels = []
    guide_labels = []
    for position, rows in selected.groupby("start_1", sort=True):
        guide_counts = rows.groupby("guide_query").size().reset_index(name="n")
        guide_counts = guide_counts.sort_values(
            ["n", "guide_query"], ascending=[False, True]
        )
        dominant_guide = str(guide_counts.iloc[0]["guide_query"])
        position_labels.append(f"`{target.star_contig}:{int(position):,}` ({len(rows):,})")
        guide_labels.append(f"`{dominant_guide}`")
    return "<br>".join(position_labels), "<br>".join(guide_labels)


def build_report(bundle: Mapping[str, object], docs_dir: Path, asset_dir: Path) -> Path:
    summary = bundle["global"]
    controlled = bundle["controlled"]
    loci = bundle["loci"].set_index("gene")
    representatives = bundle["representatives"].set_index("assigned_gene_name")
    provenance = bundle["provenance"]

    rows = []
    for target in TARGETS:
        locus = loci.loc[target.gene]
        rep = representatives.loc[target.gene]
        positions, guides = summarize_positions(controlled, target)
        rows.append(
            f"| {target.gene} | {target.gene_id} | {guides} | {positions} | {locus.annotation_region} | "
            f"{int(locus.target_concordant_reads):,} | {int(locus.removed_false_molecules):,} | "
            f"{int(locus.putative_genuine_gex_reads):,} | `{rep.contig}:{int(rep.start_1):,}` / "
            f"`{rep.cigar}` / {rep.strand} |"
        )

    sections = []
    for target in TARGETS:
        locus = loci.loc[target.gene]
        rep = representatives.loc[target.gene]
        sections.append(
            f"""### {target.gene}

Representative `{rep.truth_id}` contains guide `{rep.guide_query}` followed by
the `{rep.raw_post_guide}` stem boundary. STAR aligned the trimmed 24-mer at
`{rep.contig}:{int(rep.start_1):,}` (`{rep.cigar}`, MAPQ {int(rep.mapq)}, NM {int(rep.nm)}); the corresponding processed
full-length contaminant had no mapped primary record.

![{target.gene} contaminant and putative genuine GEX locus tracks](assets/c2_gex_grna_misalignment/{target.gene.lower()}_locus.png)
"""
        )

    report = f"""# C2 gRNA-to-GEX Misalignment Case Study

## Result

C2 contained {summary['raw_structural_contaminants']:,} full-length structural gRNA contaminants among {summary['raw_reads']:,}
raw GEX pairs ({summary['raw_structural_contaminant_cpm']:.2f} per million). The separate cell-associated pre-alignment
molecular-family estimate is {summary['prealignment_family_cpm_low']:.2f}–{summary['prealignment_family_cpm_high']:.2f} CPM; the denominators must not be
interchanged. In the controlled synthetic-54 experiment, {summary['short_guide_target_assigned_reads']:,}/{summary['input_reads']:,} known
contaminants ({format_percent(summary['guide_target_fraction_of_input'])}) were assigned to their guide's target gene. The largest
complementary examples were promoter-proximal intronic KLF2 and exonic NFKB1.

## Controlled experiment

Both arms contain the same {summary['input_reads']:,} corrected-barcode contaminants and use the
same STAR reference and options. Only the short arm truncates raw R1 to 54 nt
before identical Split_TSO preprocessing.

| Outcome | Full R1 | Synthetic 54-nt R1 |
|---|---:|---:|
| Input structural contaminants | {summary['input_reads']:,} | {summary['input_reads']:,} |
| Mapped primary reads | {summary['full_mapped_primary_reads']:,} | {summary['short_mapped_primary_reads']:,} |
| Unique MAPQ-255 reads | 0 | {summary['short_unique_reads']:,} |
| False gene-assigned reads | 0 | {summary['short_false_gene_assigned_reads']:,} |
| Intended-target assignments | 0 | {summary['short_guide_target_assigned_reads']:,} |

Target counts use the quote-normalized annotation; literal quotes in the older
controlled CSV had produced incorrect `is_guide_target=False` values.

## Representative read transformation

Both reads begin with the exact 30-nt Split_TSO
`{TSO}`. After 54-nt truncation, split-pipe removes
that TSO and sends a guide-dominated 24-mer to STAR. STAR aligns all 20 guide
bases; KLF2 also aligns one matching stem base.

The 150-nt reads resolve the CSM stem, 7-nt iBAR, partial tracr, post-tracr
bridge, and opposite-end barcode read-through. Bases 121–128 and 141–148 match
the expected reverse-complement BC1/BC2 positions; bases 129–140 are
`{R2_LINKER_2_RC}`, the reverse complement of the Parse-v3 R2 linker. This is
paired-end overlap within one library fragment, not concatenation of the R2
FASTQ record. STAR uses the same sequence for 3-prime adapter clipping.

The Query and GRCh38 letters inherit the source annotation colors: guide bases
are red, and any aligned stem-derived base is orange.

![Component-resolved raw TSO-gRNA, synthetic-54 input, and annotation-colored STAR alignment](assets/c2_gex_grna_misalignment/read_transformation.png)

| Gene | Ensembl ID | Dominant guide(s), by position | Contaminant position(s) (reads) | Parse annotation | Target-concordant reads | Removed false molecules | Putative genuine GEX reads in window | Representative STAR alignment / strand |
|---|---|---|---|---|---:|---:|---:|---|
{chr(10).join(rows)}

## Locus views

Red shows the controlled short54 contaminants; blue shows the original C2 GEX
BAM on the same GRCh38 axis. “Putative genuine GEX” denotes a primary MAPQ-255
target-`GX` read that overlaps the locus and lacks a complete residual gRNA
signature; it is not independent proof of origin. Block-aware filtering
excluded {int(loci.loc['KLF2'].excluded_span_only_locus_candidates):,} KLF2 and {int(loci.loc['NFKB1'].excluded_span_only_locus_candidates):,} NFKB1 reads whose outer spans crossed a window only
through a splice.

{chr(10).join(sections)}
## Interpretation and limitations

- The {format_percent(summary['guide_target_fraction_of_input'])} result is conditional on a known contaminant, not the whole GEX
  library.
- The contaminant-only controlled arm measures alignment potential without
  competition from genuine GEX reads. Native-library UMI collapse can merge or
  redistribute molecules, so the controlled false-molecule counts are not a
  direct correction to the expression matrix.
- The production-BAM structural classifier recognizes complete residual gRNA
  constructs. A guide-only contaminant can lack that signature, so the blue
  tracks are putative native GEX evidence, not guaranteed contamination-free
  expression.
- The original 150-cycle contaminant arm produced zero mapped reads. The
  failure appears only after shortening the same molecules and performing the
  same TSO preprocessing.
- The query-name-sorted, unindexed production BAM required one complete
  streaming pass for genuine-read recovery.

## Reproducibility

The sanitized evidence bundle in `assets/c2_gex_grna_misalignment/data/`
records UGER job `{provenance.get('job_id', 'unknown')}`, truth-panel token
`{summary['truth_panel_token']}`, provenance, sizes, and SHA-256 digests.
Read IDs are pseudonymized; cell, barcode, UMI, unused native-GEX sequence, and
absolute source paths are omitted. No STAR or full contamination scan was
repeated.

The extraction and rendering entry point is
`../scripts/c2_gex_grna_case_study.py`; the bound UGER wrapper is
`../submit/submit_c2_gex_grna_case_study.sh`. Its `publish` command validates
the internal extraction hashes and writes the scientific-minimum bundle; its
`render` command consumes only that sanitized bundle.
"""
    report_path = docs_dir / "C2_GEX_GRNA_MISALIGNMENT_CASE_STUDY.md"
    report_path.write_text(report)
    return report_path


def render_command(args: argparse.Namespace) -> None:
    data_dir = Path(args.data_dir).resolve()
    docs_dir = Path(args.docs_dir).resolve()
    asset_dir = docs_dir / "assets" / "c2_gex_grna_misalignment"
    asset_dir.mkdir(parents=True, exist_ok=True)
    bundle = load_bundle(data_dir)
    render_transformation(bundle, asset_dir)
    for target in TARGETS:
        render_locus(bundle, target, asset_dir)
    report_path = build_report(bundle, docs_dir, asset_dir)
    print(f"Case-study report rendered: {report_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract", help="Run the bounded UGER extraction")
    extract.add_argument("--config", required=True, help="Completed C2 QC config")
    extract.add_argument("--qc-root", default="", help="Override the completed QC output root")
    extract.add_argument("--target-annotation-json", default="", help="Corrected target-concordance JSON")
    extract.add_argument("--gex-bam", default="", help="Override the production annotated C2 BAM")
    extract.add_argument("--reference-fasta", default="", help="Override the indexed GRCh38 FASTA")
    extract.add_argument("--output-dir", required=True, help="Durable compact output directory")
    extract.add_argument("--scratch-dir", required=True, help="Task-specific disposable scratch")
    extract.add_argument("--threads", type=int, default=1)
    extract.add_argument("--samtools", default="")
    extract.add_argument("--pigz", default="")
    extract.set_defaults(func=extract_command)

    render = subparsers.add_parser("render", help="Render figures and Markdown from compact data")
    render.add_argument("--data-dir", required=True)
    render.add_argument("--docs-dir", default=str(PROJECT_DIR / "docs"))
    render.set_defaults(func=render_command)

    publish = subparsers.add_parser(
        "publish", help="Create a deterministic sanitized evidence bundle"
    )
    publish.add_argument("--data-dir", required=True, help="Internal extraction bundle")
    publish.add_argument(
        "--output-dir", required=True, help="Empty destination for published evidence"
    )
    publish.set_defaults(func=publish_command)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be positive")
    args.func(args)


if __name__ == "__main__":
    main()
