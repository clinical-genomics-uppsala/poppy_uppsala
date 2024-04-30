__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule bamsnap_create_pos_list:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.filter.somatic.vcf.gz",
        tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.filter.somatic.vcf.gz.tbi",
    output:
        bed=temp("bamsnap/create_pos_list/{sample}_{type}.pos.bed"),
    params:
        af=config.get("bamsnap_create_pos_list", {}).get("af", "0.05"),
        extra=config.get("bamsnap_create_pos_list", {}).get("extra", ""),
    log:
        "bamsnap/create_pos_list/{sample}_{type}.pos.bed.log",
    benchmark:
        repeat(
            "bamsnap/create_pos_list/{sample}_{type}.pos.beds.benchmark.tsv",
            config.get("bamsnap_create_pos_list", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bamsnap_create_pos_list", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap_create_pos_list", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap_create_pos_list", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap_create_pos_list", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap_create_pos_list", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap_create_pos_list", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bamsnap_create_pos_list", {}).get("container", config["default_container"])
    # localrule: True
    message:
        "{rule}: extract all positions with AF higher than {params.af} from {input.vcf}"
    script:
        "../scripts/create_pos_list.py"


rule bamsnap:
    input:
        pos_list="bamsnap/create_pos_list/{sample}_{type}.pos.bed",
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        fasta=config["reference"]["fasta"],
    output:
        results_dir=temp(directory("bamsnap/bamsnap/{sample}_{type}/")),
        index=temp("bamsnap/bamsnap/{sample}_{type}/index.html"),
        sample_list=temp("bamsnap/bamsnap/{sample}_{type}/sample_list.html"),
        variant_list=temp("bamsnap/bamsnap/{sample}_{type}/variant_list.html"),
    params:
        margin=config.get("bamsnap", {}).get("margin", "50"),
        extra=config.get("bamsnap", {}).get("extra", "-show_soft_clipped "),
    log:
        "bamsnap/bamsnap/{sample}_{type}.log",
    wildcard_constraints:
        sample="(?!HD829).*",
    benchmark:
        repeat("bamsnap/bamsnap/{sample}_{type}.benchmark.tsv", config.get("bamsnap", {}).get("benchmark_repeats", 1))
    threads: config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: create bamsnaps based on {input.pos_list} and {input.bam}"
    shell:
        "(bamsnap -bam {input.bam} -ref {input.fasta} -out {output.results_dir} -process {threads} -margin {params.margin} -bed {input.pos_list} {params.extra}) &> {log}"


rule bamsnap_hd829:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        fasta=config["reference"]["fasta"],
    output:
        results_dir=temp(directory("bamsnap/bamsnap/{sample}_{type}/")),
    log:
        "bamsnap/bamsnap/{sample}_{type}.log",
    wildcard_constraints:
        sample="(HD829).*",
    benchmark:
        repeat("bamsnap/bamsnap/{sample}_{type}.benchmark.tsv", config.get("bamsnap", {}).get("benchmark_repeats", 1))
    threads: config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: create dummy folder for {wildcards.sample} "
    shell:
        "(mkdir -p {output.results_dir}) &> {log}"
