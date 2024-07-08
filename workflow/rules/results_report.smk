__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


from hydra_genetics.utils.software_versions import get_pipeline_version


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
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz",
        vcf_tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz.tbi",
        pindel="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.artifact_annotated.filter.somatic_hard.filter.pindel.vcf.gz",
        pindel_tbi="cnv_sv/pindel_vcf/{sample}_{type}.no_tc.vep_annotated.artifact_annotated.filter.somatic_hard.filter.pindel.vcf.gz.tbi",
        pindel_bed=config["pindel_call"]["include_bed"],
        mosdepth_summary="qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.summary.txt",
        mosdepth_perbase="qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.txt",
        mosdepth_regions="qc/mosdepth_bed_coding/{sample}_{type}.regions.bed.gz",
        mosdepth_thresholds="qc/mosdepth_bed_coding/{sample}_{type}.thresholds.bed.gz",
        picard_dupl="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
        wanted_transcripts=config["results_report_xlsx"]["wanted_transcripts"],
        wait="versions/update_poppy.temp",
    output:
        xlsx="reports/xlsx/{sample}_{type}.xlsx",
    params:
        sample=lambda wildcards: wildcards.sample,
        sample_type=lambda wildcards: wildcards.type,
        sequenceid=config["sequenceid"],
        poppy_version=config["poppy_version"],
        uppsala_version=get_pipeline_version(workflow, pipeline_name="poppy_uppsala"),
        bedfile=config["reference"]["design_bed"],
        exonbed=config["reference"]["exon_bed"],
        pindelbed=config["pindel_call"]["include_bed"],
        ref=config["reference"]["fasta"],
        artifact=config["reference"]["artifacts"],
        background=config["reference"]["background"],
        artifact_pindel=config["reference"]["artifacts_pindel"],
        thresholds=config["mosdepth_bed"]["thresholds"],
        extra=config.get("results_report", {}).get("extra", ""),
    log:
        "reports/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat(
            "reports/xlsx/{sample}_{type}.xlsx.benchmark.tsv", config.get("results_report", {}).get("benchmark_repeats", 1)
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
        "{rule}: summerize results into {output.xlsx}."
    # localrule: True
    script:
        "../scripts/results_report_xlsx.py"
