#!/bin/python3

from pysam import VariantFile

vcf=VariantFile(snakemake.input.vcf)

with open(snakemake.output.bed, "w+") as output:
    for record in vcf.fetch():
        if record.filter.keys() == ["PASS"] and record.info["AF"][0] >= float(snakemake.params.af):
            output.write("\t".join([record.contig, str(record.pos), str(record.pos)])+"\n")