#!/usr/bin/env python3
"""Run beside the four input files: python cre_te_simple_simplified.py"""

import csv
import gzip
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    from pybedtools import BedTool, helpers
except ImportError:
    sys.exit("Install pybedtools first: pip install pybedtools")


# Change these four names only if your input filenames are different.
HERE = path(__file__).parent
CRE_BED = HERE/"bed/HEK_Neuron_THP1.CRE.coord.bed.gz"
CRE_INFO = HERE/"log/HEK_Neuron_THP1.CRE.info.p.e.se.tsv"
COUNTS = HERE/"counts/HEK_Neuron_THP1.all.counts.tsv"
RMSK = "/group/carninci/Shakiba.hamedi/CRE_all/analysis/CRE_all/rmsk.txt.gz"
OUTDIR = HERE/"cre_te_results"

PREFIXES = (
    ("HEK_", "HEK_RPI"), ("iPSC_", "iPSC"), ("NSC_", "NSC"),
    ("Neuron_", "Neuron"), ("THP1_", "THP1"),
)
CELL_TYPES = tuple(name for _, name in PREFIXES)
BED_HEADER = "chrom start end CREID cre_score cre_strand summit_start summit_end".split()
EXTRA_HEADER = "cell_type cell_type_total_count nonzero_sample_count info_record_found counts_record_found".split()
TE_HEADER = "TE_chrom TE_start TE_end TE_name TE_class TE_family TE_orientation TE_id bp_overlap pct_CRE_overlap pct_TE_overlap overlap_status".split()
TE_CLASSES = {"DNA", "LINE", "LTR", "SINE", "RC", "Retroposon"}


def cell_type(sample):
    return next((name for prefix, name in PREFIXES if sample.startswith(prefix)), None)


def percent(n, total):
    return f"{100 * n / total:.4f}" if total else "NA"


def read_info(path):
    with path.open(newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = ["cre_orientation" if x == "orientation" else x for x in next(rows)[1:]]
        return header, {row[0]: row[1:] for row in rows if row and row[0]}


def read_counts(path):
    counts = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    with path.open(newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        types = [cell_type(sample) for sample in next(rows)[1:]]
        for row in rows:
            if not row or not row[0]:
                continue
            for kind, value in zip(types, row[1:]):
                if kind:
                    value = int(value or 0)
                    counts[row[0]][kind][0] += value
                    counts[row[0]][kind][1] += value > 0
    return counts


def master_rows(bed, info, counts, info_header):
    cre_id = bed[3]
    found = "yes" if cre_id in info else "no"
    base = bed[:8] + info.get(cre_id, [""] * len(info_header))
    if cre_id not in counts:
        return [base + ["NA", ".", ".", found, "no"]]
    rows = []
    for kind in CELL_TYPES:
        total, nonzero = counts[cre_id].get(kind, [0, 0])
        if total:
            rows.append(base + [kind, str(total), str(nonzero), found, "yes"])
    return rows or [base + ["none", "0", "0", found, "yes"]]


def write_master(info_header, info, counts, master_path, master_bed):
    header = BED_HEADER + info_header + EXTRA_HEADER
    with gzip.open(CRE_BED, "rt", newline="") as source, \
            master_path.open("w", newline="") as master_file, \
            master_bed.open("w", newline="") as bed_file:
        master_writer = csv.writer(master_file, delimiter="\t", lineterminator="\n")
        bed_writer = csv.writer(bed_file, delimiter="\t", lineterminator="\n")
        master_writer.writerow(header)

        for bed in csv.reader(source, delimiter="\t"):
            if len(bed) < 12:
                continue
            for row in master_rows(bed, info, counts, info_header):
                master_writer.writerow(row)
                bed_writer.writerow(row)
    return header


def write_repeatmasker_bed(path, te_bed):
    with gzip.open(path, "rt", newline="") as source, te_bed.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        for row in csv.reader(source, delimiter="\t"):
            if len(row) < 17:
                continue
            writer.writerow([
                row[5], row[6], row[7], row[10], row[11], row[12],
                "-" if row[9] == "C" else row[9], row[16],
            ])


def write_te_overlaps(master_bed, te_bed, master_header, output_path):
    overlaps = BedTool(str(master_bed)).sort().intersect(
        BedTool(str(te_bed)).sort(), wao=True, sorted=True
    )
    width = len(master_header)

    with Path(overlaps.fn).open(newline="") as source, \
            output_path.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(master_header + TE_HEADER)
        for row in csv.reader(source, delimiter="\t"):
            master, te, bp = row[:width], row[width:-1], int(row[-1])
            if not bp or te[0] == ".":
                continue
            te[4] = te[4].rstrip("?")  # LINE? -> LINE, SINE? -> SINE
            if te[4] not in TE_CLASSES:
                continue
            te_length = int(te[2]) - int(te[1])
            cre_length = int(master[2]) - int(master[1])
            writer.writerow(master + te + [
                str(bp), percent(bp, cre_length), percent(bp, te_length), "overlap"
            ])


def main():
    if not shutil.which("bedtools"):
        sys.exit("Install bedtools first: conda install -c bioconda bedtools")

    OUTDIR.mkdir(exist_ok=True)
    info_header, info = read_info(CRE_INFO)
    counts = read_counts(COUNTS)
    master_path = OUTDIR / "CRE_master.tsv"
    overlap_path = OUTDIR / "CRE_TE_overlaps.tsv"

    with tempfile.TemporaryDirectory(dir=OUTDIR, prefix=".tmp_") as temp_dir:
        temp_dir = Path(temp_dir)
        helpers.set_tempdir(str(temp_dir))
        master_bed = temp_dir / "master.bed"
        te_bed = temp_dir / "repeatmasker.bed"

        header = write_master(info_header, info, counts, master_path, master_bed)
        write_repeatmasker_bed(RMSK, te_bed)
        write_te_overlaps(master_bed, te_bed, header, overlap_path)

    print("Done:", master_path, overlap_path, sep="\n")


if __name__ == "__main__":
    main()
