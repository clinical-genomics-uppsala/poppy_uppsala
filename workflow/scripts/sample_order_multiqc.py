#!/bin/python3

import sys
import csv
import re

# process input as pairs, [sample, fastq1] from units.tsv
sample_order_duplicates = []
for sample, fastq_path in snakemake.params.filelist:
    fastq = fastq_path.split("/")[-1]
    s_pattern = re.compile("_S([0-9]+)_")
    s_index = int(s_pattern.search(fastq).group(1))
    sample_order_duplicates.append([sample, s_index])


# Remove duplicates from sample_order_duplicates and order based on s_index in fastq1 filename
sample_order = [list(x) for x in set(tuple(x) for x in sample_order_duplicates)]
sample_order.sort(key=lambda x: int(x[1]))


with open(snakemake.output.replacement, "w+") as tsv:
    tsv_writer = csv.writer(tsv, delimiter="\t")
    i = 1
    for sample in sample_order:
        tsv_writer.writerow([sample[0], "sample_" + str(f"{i:03}")])

        i += 1

with open(snakemake.output.order, "w+") as tsv:
    tsv_writer = csv.writer(tsv, delimiter="\t")
    tsv_writer.writerow(["Sample Order", "Sample Name"])
    i = 1
    for sample in sample_order:
        tsv_writer.writerow(["sample_" + str(f"{i:03}"), sample[0]])
        i += 1
