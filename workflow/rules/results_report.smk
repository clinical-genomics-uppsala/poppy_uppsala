__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule results_report_xlsx:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.filter.somatic.vcf",
        pindel="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.filtered.pindel.vcf.gz",
    output:
        xlsx="results_report/xlsx/{sample}_{type}.xlsx",
    params:
        sequenceid=config["sequenceid"],
        poppy_version="test",
        uppsala_version="test",
        ref=config["references"]["fasta"],
        pindel_bed=config["pindel_call"]["include_bed"],
        extra=config.get("results_report", {}).get("extra", ""),
    log:
        "results_report/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "results_report/xlsx/{sample}_{type}.xlsx.benchmark.tsv",
            config.get("results_report", {}).get("benchmark_repeats", 1)
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
