__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2026, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnv_large_calls_table:
    input:
        vcf="cnv_sv/svdb_query/{sample}_{type}.{tc_method}.svdb_query.annotate_cnv.cnv_genes.vcf.gz",
        cytobands=config["merge_cnv_json"]["cytobands"],
    output:
        tsv=temp("reports/cnv_html_report/{sample}_{type}.{tc_method}.large_cnvs.tsv"),
    params:
        min_length=config.get("cnv_large_calls_table", {}).get("min_length", 100000),
        max_normal_af=config.get("cnv_large_calls_table", {}).get("max_normal_af", 0.15),
    log:
        "reports/cnv_html_report/{sample}_{type}.{tc_method}.large_cnvs.tsv.log",
    benchmark:
        repeat(
            "reports/cnv_html_report/{sample}_{type}.{tc_method}.large_cnvs.tsv.benchmark.tsv",
            config.get("cnv_large_calls_table", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnv_large_calls_table", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnv_large_calls_table", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_large_calls_table", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_large_calls_table", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnv_large_calls_table", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnv_large_calls_table", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_large_calls_table", {}).get("container", config["default_container"])
    message:
        "{rule}: list CNVs larger than {params.min_length} bp in {output.tsv}"
    script:
        "../scripts/cnv_large_calls_table.py"
