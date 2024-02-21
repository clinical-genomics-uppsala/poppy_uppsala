__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule results_report_bedtools_intersect:
    input:
        left="qc/mosdepth_bed_coding/{sample}_{type}.per-base.bed.gz",
        coverage_csi="qc/mosdepth_bed_coding/{sample}_{type}.per-base.bed.gz.csi",
        right=config["reference"]["exon_bed"],
    output:
        results=temp("qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.txt"),
    params:
        extra=config.get("results_report_bedtools_intersect", {}).get("extra", ""),
    log:
        "qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.log",
    benchmark:
        repeat(
            "qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.benchmark.tsv",
            config.get("results_report_bedtools_intersect", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("results_report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("results_report_bedtools_intersect", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("results_report_bedtools_intersect", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("results_report_bedtools_intersect", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("results_report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("results_report_bedtools_intersect", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("results_report_bedtools_intersect", {}).get("container", config["default_container"])
    message:
        "{rule}: export low cov regions from {input.left} based on {input.right}"
    wrapper:
        "v1.32.0/bio/bedtools/intersect"


rule results_report_xlsx:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.filter.somatic.vcf.gz",
        vcf_tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.filter.somatic.vcf.gz.tbi",
        pindel="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.vcf.gz",
        pindel_tbi="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.vcf.gz.tbi",
        pindel_bed=config["pindel_call"]["include_bed"],
        mosdepth_perbase="qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.txt",
        mosdepth_regions="qc/mosdepth_bed_coding/{sample}_{type}.regions.bed.gz",
        # wanted_transcripts=config["xlsx_report"]["wanted_transcripts"],
    output:
        xlsx="results_report/xlsx/{sample}_{type}.xlsx",
    params:
        sample=lambda wildcards: wildcards.sample,
        sequenceid=config["sequenceid"],
        poppy_version="test",
        uppsala_version="test",
        bedfile=config["reference"]["design_bed"],
        exonbed=config["reference"]["exon_bed"],
        pindelbed=config["pindel_call"]["include_bed"],
        ref=config["reference"]["fasta"],
        artifact=config["reference"]["artifacts"],
        extra=config.get("results_report", {}).get("extra", ""),
    log:
        "results_report/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "results_report/xlsx/{sample}_{type}.xlsx.benchmark.tsv", config.get("results_report", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("results_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("results_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("results_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("results_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("results_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("results_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("results_report", {}).get("container", config["default_container"])
    message:
        "{rule}: summerize results into {output.xlsx}"
    script:
        "../scripts/results_report_xlsx.py"
