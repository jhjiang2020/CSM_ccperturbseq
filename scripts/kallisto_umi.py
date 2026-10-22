"""Streaming native BUS counting with Parse 1.6.4 UMI-collapse semantics.

No barcode correction or BAM projection occurs here. Corrected barcode-index
identities remain distinct until after molecule counting.
"""

from __future__ import annotations

import gzip
import itertools
import sqlite3
import struct
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


BASES = "ACGT"
BUS_RECORD = struct.Struct("<QQiIII")
UMI_POLICY = "parse_1.6.4_strict_abundance_hamming1_v1"
UMI_POLICIES = {
    "parse": UMI_POLICY,
    "directional": "directional_hamming1_2x_minus1_bustools_traversal_v1",
    "none": "exact_umi_only_v1",
}


def umi_policy(flavor):
    if flavor not in UMI_POLICIES:
        raise ValueError("umi_correction must be parse, directional, or none")
    return UMI_POLICIES[flavor]


def dna_encode(sequence):
    result = 0
    for base in sequence:
        if base not in BASES:
            raise ValueError(f"Non-ACGT barcode/UMI: {sequence!r}")
        result = (result << 2) | BASES.index(base)
    return result


def dna_decode(value, length):
    return "".join(BASES[(value >> (2 * i)) & 3] for i in reversed(range(length)))


def fastq_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline().rstrip("\r\n")
            if not header:
                return
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline().rstrip("\r\n")
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header, sequence, quality


def prepare_fastqs(source, cdna_path, tag_path):
    """Reformat Parse survivors; the barcode is the *uncollapsed* BC sequence.

    Return barcode lookup, barcode length, UMI length, and read accounting.
    N-containing UMIs cannot be represented losslessly in native BUS and are
    excluded explicitly, rather than silently converted by kallisto.
    """
    lookup = {}
    lengths = None
    stats = Counter(input_reads=0, submitted_reads=0, invalid_umi_reads=0, empty_cdna_reads=0)
    with open(cdna_path, "w") as cdna, open(tag_path, "w") as tags:
        for header, sequence, quality in fastq_records(source):
            stats["input_reads"] += 1
            parts = header[1:].split()[0].split("__")
            if len(parts) < 5:
                raise ValueError("Expected Parse barcode_head FASTQ headers")
            wells, primer, bcis, barcode_parts, umi = parts[:5]
            barcode = barcode_parts.replace("_", "")
            if len(wells.split("_")) != 3 or len(bcis.split("_")) != 3:
                raise ValueError(f"Malformed Parse barcode identity: {header}")
            dna_encode(barcode)
            current_lengths = (len(barcode), len(umi))
            if not all(0 < value <= 32 for value in current_lengths):
                raise ValueError("BUS requires barcode and UMI lengths from 1 to 32")
            if lengths is not None and lengths != current_lengths:
                raise ValueError("Inconsistent barcode/UMI lengths in Parse survivors")
            lengths = current_lengths
            identity = (wells, bcis, primer, barcode_parts)
            if barcode in lookup and lookup[barcode] != identity:
                raise ValueError("One corrected barcode maps to conflicting Parse identities")
            lookup[barcode] = identity
            if any(base not in BASES for base in umi):
                stats["invalid_umi_reads"] += 1
                continue
            if not sequence:
                stats["empty_cdna_reads"] += 1
                continue
            name = f"@csm{stats['submitted_reads']}"
            cdna.write(f"{name}\n{sequence}\n+\n{quality}\n")
            tag = barcode + umi
            tags.write(f"{name}\n{tag}\n+\n{'I' * len(tag)}\n")
            stats["submitted_reads"] += 1
    if lengths is None or not stats["submitted_reads"]:
        raise ValueError("No BUS-representable transcriptomic reads after preprocessing")
    return lookup, *lengths, dict(stats)


def bus_records(path, barcode_length, umi_length):
    """Read BUS v1 records with strict validation (including short final records)."""
    with open(path, "rb") as handle:
        header = handle.read(20)
        if len(header) != 20:
            raise ValueError("Truncated BUS header")
        magic, version, bc_len, umi_len, text_length = struct.unpack("<4sIIII", header)
        if magic != b"BUS\x00" or version != 1:
            raise ValueError("Expected uncompressed BUS v1")
        if (bc_len, umi_len) != (barcode_length, umi_length):
            raise ValueError("BUS barcode/UMI lengths disagree with reformatted input")
        if len(handle.read(text_length)) != text_length:
            raise ValueError("Truncated BUS header text")
        barcodes = {}
        while True:
            record = handle.read(BUS_RECORD.size)
            if not record:
                return
            if len(record) != BUS_RECORD.size:
                raise ValueError("Truncated BUS record")
            bc, umi, ec, count, flags, padding = BUS_RECORD.unpack(record)
            if ec < 0 or count < 1 or flags != 0 or padding != 0:
                raise ValueError("Invalid BUS record or unsupported flags")
            if bc >> (2 * bc_len) or umi >> (2 * umi_len):
                raise ValueError("BUS barcode/UMI exceeds its declared length")
            barcode = barcodes.get(bc)
            if barcode is None:
                barcode = barcodes[bc] = dna_decode(bc, bc_len)
            yield barcode, dna_decode(umi, umi_len), ec, count


def read_equivalence_classes(path, n_transcripts):
    classes = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 2:
                raise ValueError("Malformed kallisto matrix.ec")
            ec = int(fields[0])
            targets = tuple(sorted(set(map(int, fields[1].split(",")))))
            if ec in classes or not targets or min(targets) < 0 or max(targets) >= n_transcripts:
                raise ValueError("Invalid or duplicated equivalence class")
            classes[ec] = targets
    return classes


def spool_bus(bus_path, database, lookup, classes, bc_len, umi_len):
    """Disk-backed aggregation; no cell/UMI group is split across workers/lanes."""
    connection = sqlite3.connect(database)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute("CREATE TABLE reads (barcode TEXT, umi TEXT, ec INTEGER, n INTEGER)")
    batch = []
    read_count = 0
    for barcode, umi, ec, count in bus_records(bus_path, bc_len, umi_len):
        if barcode not in lookup or ec not in classes:
            raise ValueError("BUS barcode or equivalence class is absent from input/reference")
        batch.append((barcode, umi, ec, count))
        read_count += count
        if len(batch) == 100_000:
            connection.executemany("INSERT INTO reads VALUES (?,?,?,?)", batch)
            batch.clear()
    connection.executemany("INSERT INTO reads VALUES (?,?,?,?)", batch)
    connection.execute("CREATE INDEX reads_order ON reads(barcode,umi,ec)")
    connection.commit()
    return connection, read_count


def single_mutations(sequence):
    for index, old_base in enumerate(sequence):
        for base in BASES:
            if base != old_base:
                yield sequence[:index] + base + sequence[index + 1:]


def directional_umis(read_counts):
    """Fixed-count directional graph, with bustools' root/ownership ordering.

    Correction is scoped to one barcode and one candidate-gene set by callers.
    Edges use original read support, including bidirectional singleton edges.
    Roots claim reachable descendants transitively in lexical UMI order.
    """
    ordered = sorted(read_counts)
    if len({len(umi) for umi in ordered}) > 1:
        raise ValueError("UMIs in one group have different lengths")
    if len(ordered) < 2:
        return [(umi, read_counts[umi], (umi,)) for umi in ordered]
    outgoing = {umi: set() for umi in ordered}
    incoming = {umi: set() for umi in ordered}
    for umi in ordered:
        for neighbor in single_mutations(umi):
            if neighbor not in read_counts:
                continue
            for parent, child in ((umi, neighbor), (neighbor, umi)):
                if read_counts[parent] >= 2 * read_counts[child] - 1:
                    outgoing[parent].add(child)
                    incoming[child].add(parent)
    remaining = set(ordered)
    result = []
    for first in ordered:
        if first not in remaining:
            continue
        component, stack = set(), [first]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend((outgoing[node] | incoming[node]) - component)
        pending = component & remaining
        while pending:
            roots = sorted(node for node in pending if not incoming[node]) or [min(pending)]
            for root in roots:
                members, stack = set(), [root]
                while stack:
                    node = stack.pop()
                    if node not in pending or node in members:
                        continue
                    members.add(node)
                    stack.extend(outgoing[node] - members)
                pending.difference_update(members)
                remaining.difference_update(members)
                result.append((root, sum(read_counts[node] for node in members), tuple(sorted(members))))
    return result


def collapse_umis(read_counts, *, flavor="parse", kernel=None, mutation_function=single_mutations):
    """Parse's strict-abundance algorithm, retaining merged-family membership.

    Counts are sorted once. The representative's comparison count is captured
    before visiting its neighbors, exactly as in Parse (no directional ratio,
    no transitive graph closure, and no merging of equal abundance neighbors).
    When supplied, the installed Parse kernel verifies the result directly.
    """
    umi_policy(flavor)
    if flavor == "none":
        return [(umi, read_counts[umi], (umi,)) for umi in sorted(read_counts)]
    if flavor == "directional":
        return directional_umis(read_counts)
    ordered = sorted(read_counts, key=lambda umi: -read_counts[umi])
    neighbors = defaultdict(list)
    if len(ordered) > 50:
        present = set(ordered)
        for umi in ordered:
            for neighbor in mutation_function(umi):
                if neighbor in present:
                    neighbors[neighbor].append(umi)
    else:
        for first, second in itertools.combinations(ordered, 2):
            if len(first) != len(second):
                raise ValueError("UMIs in one group have different lengths")
            if sum(a != b for a, b in zip(first, second)) <= 1:
                neighbors[first].append(second)
                neighbors[second].append(first)
    counts = {umi: read_counts[umi] for umi in ordered}
    members = {umi: [umi] for umi in ordered}
    for umi in ordered:
        if umi not in counts:
            continue
        original_count = counts[umi]
        for neighbor in neighbors[umi]:
            if neighbor in counts and counts[neighbor] < original_count:
                counts[umi] += counts.pop(neighbor)
                members[umi].extend(members.pop(neighbor))
    if kernel is not None:
        expected = kernel(ordered, [read_counts[umi] for umi in ordered])
        if counts != expected:
            raise RuntimeError("Custom UMI families disagree with installed Parse collapse kernel")
    return [(umi, count, tuple(members[umi])) for umi, count in counts.items()]


@dataclass
class Evidence:
    umi: str
    genes: frozenset[str]
    transcripts: frozenset[int]
    reads: int


def resolve_exact_umi(evidence):
    """Partition UMI collisions while intersecting compatible EC evidence.

    Singleton gene anchors are retained independently. An ambiguous record
    supporting existing anchors adds read support without inventing a new
    molecule. Disjoint gene evidence represents separate molecules.
    """
    groups = []
    for item in sorted(evidence, key=lambda item: (len(item.genes), tuple(sorted(item.genes)), tuple(sorted(item.transcripts)))):
        matches = [group for group in groups if group.genes & item.genes]
        if len(matches) > 1:
            # No principled read allocation among independent anchored molecules.
            # Keep existing support; do not count this ambiguous evidence again.
            continue
        if matches:
            group = matches[0]
            group.genes &= item.genes
            group.transcripts &= item.transcripts
            group.reads += item.reads
        else:
            groups.append(Evidence(item.umi, item.genes, item.transcripts, item.reads))
    return groups


def count_barcode(rows, classes, transcript_genes, *, fractional=False,
                  min_reads=1, kernel=None, mutation_function=single_mutations,
                  umi_correction="parse", class_evidence=None):
    """Yield (representative UMI, candidate genes, read support, transcript set)."""
    by_candidates = defaultdict(dict)
    for umi, records in itertools.groupby(rows, key=lambda row: row[0]):
        evidence = []
        for _, ec, reads in records:
            if class_evidence is None:
                targets = frozenset(classes[ec])
                genes = frozenset(transcript_genes[target] for target in targets)
            else:
                genes, targets = class_evidence[ec]
            if not fractional and len(genes) != 1:
                continue
            evidence.append(Evidence(umi, genes, targets, reads))
        resolved = resolve_exact_umi(evidence) if fractional else evidence
        for item in resolved:
            group = by_candidates[tuple(sorted(item.genes))]
            if umi not in group:
                group[umi] = Evidence(umi, item.genes, item.transcripts, item.reads)
            else:
                group[umi].reads += item.reads
                group[umi].transcripts &= item.transcripts
    for genes, items in sorted(by_candidates.items()):
        read_counts = {umi: item.reads for umi, item in items.items()}
        for umi, reads, members in collapse_umis(
            read_counts, flavor=umi_correction, kernel=kernel, mutation_function=mutation_function
        ):
            if reads < min_reads:
                continue
            targets = set(items[members[0]].transcripts)
            for member in members[1:]:
                targets.intersection_update(items[member].transcripts)
            yield umi, genes, reads, tuple(sorted(targets))
