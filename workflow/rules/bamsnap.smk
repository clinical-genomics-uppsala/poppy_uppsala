__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2024, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule bamsnap_create_pos_list:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz",
        tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.background_annotated.filter.somatic_hard.filter.somatic.vcf.gz.tbi",
    output:
        bed=temp("bamsnap/create_pos_list/{sample}_{type}.pos.bed"),
    params:
        af=config.get("bamsnap_create_pos_list", {}).get("af", "0.05"),
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


rule bamsnap_samtools_view_dedup:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
    output:
        bam=temp("bamsnap/bamsnap_samtools_view_dedup/{sample}_{type}.bam"),
    params:
        extra=config.get("bamsnap_samtools_view_dedup", {}).get("extra", ""),
    log:
        "bamsnap/bamsnap_samtools_view_dedup/{sample}_{type}.bam.log",
    benchmark:
        repeat(
            "bamsnap/bamsnap_samtools_view_dedup/{sample}_{type}.bam.benchmark.tsv",
            config.get("bamsnap_samtools_view_dedup", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bamsnap_samtools_view_dedup", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap_samtools_view_dedup", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap_samtools_view_dedup", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap_samtools_view_dedup", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap_samtools_view_dedup", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap_samtools_view_dedup", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bamsnap_samtools_view_dedup", {}).get("container", config["default_container"])
    message:
        "{rule}: create bam {output} with deduplicated reads reads from {input.bam}"
    shell:
        "(samtools view -@ {threads} -F 1024 {params.extra} -b {input.bam} > {output.bam}) &> {log}"


rule bamsnap_downsample_bam:
    input:
        bam="bamsnap/bamsnap_samtools_view_dedup/{sample}_{type}.bam",
    output:
        bam=temp("bamsnap/bamsnap_downsample_bam/{sample}_{type}.bam"),
    params:
        extra_count=config.get("bamsnap_downsample_bam", {}).get("extra_count", ""),
        extra_subsample=config.get("bamsnap_downsample_bam", {}).get("extra_subsample", ""),
        filter_reads=config.get("bamsnap_downsample_bam", {}).get("filter_reads", "-F2060"),
        max_reads=config.get("bamsnap_downsample_bam", {}).get("max_reads", 2500000),
        float_precision=config.get("bamsnap_downsample_bam", {}).get("float_precision", 3),
    log:
        "bamsnap/bamsnap_downsample_bam/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "bamsnap/bamsnap_downsample_bam/{sample}_{type}.output.benchmark.tsv",
            config.get("bamsnap_downsample_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bamsnap_downsample_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bamsnap_downsample_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bamsnap_downsample_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bamsnap_downsample_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bamsnap_downsample_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bamsnap_downsample_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bamsnap_downsample_bam", {}).get("container", config["default_container"])
    message:
        "{rule}: Output only a proportion of the alignments in {input} based on a maximum number of reads equal to {params.max_reads}"
    shell:
        "nb_reads=$(samtools view -c {params.filter_reads} {params.extra_count} {input.bam}) &> {log} && "
        "frac_reads=$( bc -l <<< \"scale={params.float_precision}; ({params.max_reads}-1)/${{nb_reads}}\" ) &>> {log} "
        "&& "
        "if (( $( bc -l <<< \"$frac_reads < 1\" ) )); then "
        "echo \"File has more than {params.max_reads} reads, downsampling to ca. {params.max_reads} reads (fraction: "
        "$frac_reads)\" &>> {log} && "
        "(samtools view -@ {threads} --subsample $frac_reads {params.extra_subsample} -b {input.bam} > {output.bam}) "
        "&>> {log}; "
        "else echo \"File has fewer than {params.max_reads} reads, no subsampling was done.\" &>> {log} && "
        "(samtools view -@ {threads} -b {input.bam} > {output.bam}) &>> {log}; fi"


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
    container:
        config.get("bamsnap", {}).get("container", config["default_container"])
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
    # localrule: True
    container:
        config.get("bamsnap_hd829", {}).get("container", config["default_container"])
    message:
        "{rule}: create dummy folder for {wildcards.sample} "
    shell:
        "(mkdir -p {output.results_dir}) &> {log}"

