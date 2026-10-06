__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2026, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import csv
import sys

from pysam import VariantFile

HEADER = ["Chromosome", "Start", "End", "Cytoband", "Type", "Copy number", "Length", "Caller"]


def read_cytobands(filename):
    """Read a UCSC cytoBand file (0-based, half-open) into {chrom: [(start, end, band), ...]}."""
    cytobands = {}
    with open(filename) as f:
        for line in f:
            chrom, start, end, band = line.rstrip("\n").split("\t")[:4]
            cytobands.setdefault(chrom, []).append((int(start), int(end), band))
    for bands in cytobands.values():
        bands.sort()
    return cytobands


def get_cytoband(cytobands, chrom, start, end):
    """Return the cytoband range, e.g. '1p36.33-p36.31', overlapping the 1-based closed interval start-end."""
    bands = [band for band_start, band_end, band in cytobands.get(chrom, []) if band_start < end and band_end >= start]
    if not bands:
        return ""
    chrom_name = chrom.removeprefix("chr")
    if len(bands) == 1:
        return f"{chrom_name}{bands[0]}"
    return f"{chrom_name}{bands[0]}-{bands[-1]}"


def first_value(value):
    if isinstance(value, (list, tuple)):
        return value[0]
    return value


def get_large_cnvs(vcf_filename, cytobands, min_length, max_normal_af):
    """Return table rows for gains and losses longer than min_length that are not common in the normal panel."""
    rows = []
    for record in VariantFile(vcf_filename):
        svtype = record.info.get("SVTYPE")
        length = abs(first_value(record.info.get("SVLEN")))
        normal_af = first_value(record.info.get("Normal_AF", 0))
        if svtype == "COPY_NORMAL" or length <= min_length or normal_af > max_normal_af:
            continue
        copy_number = first_value(record.info.get("CORR_CN"))
        rows.append(
            {
                "Chromosome": record.chrom,
                "Start": record.pos,
                "End": record.stop,
                "Cytoband": get_cytoband(cytobands, record.chrom, record.pos, record.stop),
                "Type": svtype,
                "Copy number": "" if copy_number is None else round(copy_number, 2),
                "Length": length,
                "Caller": record.info.get("CALLER"),
            }
        )
    return rows


def write_table(rows, filename):
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w")
    cytobands = read_cytobands(snakemake.input.cytobands)
    rows = get_large_cnvs(snakemake.input.vcf, cytobands, snakemake.params.min_length, snakemake.params.max_normal_af)
    write_table(rows, snakemake.output.tsv)


if __name__ == "__main__":
    main()
