#!/usr/bin/env python3
"""Build a checksum-pinned kallisto mature+nascent reference using kb ref.

Run on a compute node. This creates a new directory and never changes Parse.
"""

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

try:
    from scripts.gex_alignment import SCHEMA, atomic_json, sha256_file
except ModuleNotFoundError:
    from gex_alignment import SCHEMA, atomic_json, sha256_file


def feature_annotations(gtf):
    annotations = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed GTF")
            attributes = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
            gene = attributes.get("gene_id")
            if gene:
                name = attributes.get("gene_name", gene)
                # kb/ngs-tools use_version=True adds the *explicit* GTF version.
                # Resolve that alias to gene_id; never strip arbitrary suffixes.
                for alias in (gene, f"{gene}.{attributes['gene_version']}" if attributes.get("gene_version") else gene):
                    value = (gene, name)
                    if alias in annotations and annotations[alias] != value:
                        raise ValueError(f"Conflicting GTF gene alias: {alias}")
                    annotations[alias] = value
    if not annotations:
        raise ValueError("GTF has no annotated genes")
    return annotations


def finalize_reference(args, temporary, command):
    """Validate kb runtime assets and publish their annotated reference contract."""
    temporary = Path(temporary)
    output = Path(args.output).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Reference destination already exists: {output}")
    fasta, gtf = Path(args.fasta).resolve(), Path(args.gtf).resolve()
    # Verify the configured runtime can read the built index.
    subprocess.run([args.kallisto, "inspect", "-t", str(args.threads), str(temporary / "index.idx")], check=True,
                   stdout=subprocess.DEVNULL)
    nascent_path = temporary / "nascent.txt"
    nascent = set(nascent_path.read_text().splitlines()) if nascent_path.exists() else set()
    names = feature_annotations(gtf)
    seen = set()
    canonical_t2g = []
    with (temporary / "features.tsv").open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_name", "genome", "feature_type"])
        with (temporary / "t2g.txt").open() as t2g:
            for line in t2g:
                fields = line.rstrip().split("\t")
                transcript, kb_gene = fields[:2]
                gene, name = names.get(kb_gene, (None, None))
                if transcript in seen or gene is None:
                    raise ValueError(f"Duplicated target or unannotated gene: {transcript}, {kb_gene}")
                seen.add(transcript)
                fields[1] = gene
                canonical_t2g.append("\t".join(fields))
                writer.writerow([transcript, gene, name, args.genome,
                                 "nascent" if transcript in nascent else "transcript"])
    (temporary / "t2g.txt").write_text("\n".join(canonical_t2g) + "\n")
    files = ["index.idx", "t2g.txt", "features.tsv"]
    files += [name for name in ("cdna.txt", "nascent.txt") if (temporary / name).exists()]
    version = subprocess.check_output([args.kallisto, "version"], text=True).strip()
    manifest = {
        "schema_version": SCHEMA, "complete": True, "workflow": args.workflow,
        "kallisto_version": version, "command": command,
        "gene_id_policy": "gtf_gene_id_with_explicit_version_aliases",
        "sources": {"fasta_sha256": sha256_file(fasta), "gtf_sha256": sha256_file(gtf)},
        "artifacts": {name: sha256_file(temporary / name) for name in files},
    }
    atomic_json(temporary / "reference.json", manifest)
    # Only runtime assets and the build log become durable reference files.
    output.mkdir()
    for name in files + ["reference.json", "build.log"]:
        shutil.move(str(temporary / name), output / name)
    return output


def build(args):
    output = Path(args.output).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Reference destination already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fasta, gtf = Path(args.fasta).resolve(), Path(args.gtf).resolve()
    for path in (fasta, gtf):
        if not path.is_file():
            raise FileNotFoundError(path)
    with tempfile.TemporaryDirectory(prefix=".kallisto-reference-", dir=output.parent) as tmp:
        temporary = Path(tmp)
        command = [args.kb, "ref", "--workflow", args.workflow,
                   "--kallisto", args.kallisto, "-i", str(temporary / "index.idx"),
                   "-g", str(temporary / "t2g.txt"), "-f1", str(temporary / "cdna.fa"),
                   "-t", str(args.threads)]
        if args.workflow == "nac":
            command += ["-f2", str(temporary / "nascent.fa"),
                        "-c1", str(temporary / "cdna.txt"), "-c2", str(temporary / "nascent.txt")]
        command += ["--tmp", str(temporary / "kb_tmp"), str(fasta), str(gtf)]
        with (temporary / "build.log").open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            failed_log = output.with_name(output.name + ".failed.log")
            shutil.copy2(temporary / "build.log", failed_log)
            raise RuntimeError(f"kb ref failed with exit code {result.returncode}; see {failed_log}")
        finalize_reference(args, temporary, command)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome", required=True, help="Exact Parse genome label, e.g. hg38")
    parser.add_argument("--output", required=True)
    parser.add_argument("--kb", default="kb")
    parser.add_argument("--kallisto", default="kallisto")
    parser.add_argument("--workflow", choices=("nac", "standard"), default="nac")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    print(build(args))


if __name__ == "__main__":
    main()
