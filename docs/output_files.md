# Output files and delivered files

Many temporary files, only a subset of them is forwarded to the geneticists.
These final output files are listed and shortly described in the table below.
`{sample}` is replaced with the sample name and `{type}` is replaced with the sample type (T/N/R, where T=Tumor, N=Normal, R=RNA). 
`{analysisID}` is replaced with the ID of the current analysis.

The **important files** that are delivered to the geneticists are shown in **bold**.

## Structure of the folder with the delivered files

The following files are delivered for the analysis in `results/`:

```bash
results/
  <sampleID>_<type>_<analysisID>/
    data/
      <sampleID>_<type>_<analysisID>.bam
      <sampleID>_<type>_<analysisID>.bam.bai
      <sampleID>_<type>_<analysisID>.pathology.svdb_query.vcf.gz
      <sampleID>_<type>_<analysisID>.snv-indels.vcf.gz
    reports/
      <sampleID>_<type>_<analysisID>_bamsnap/
      <sampleID>_<type>_<analysisID>_bamsnap
      
  multiqc_<analysisID>.html
```

That is, there is a separate folder for each sample.

## Data output files
These files are used by the geneticists whenever some extra checks are needed under the analysis of the variants.
<br />

| **File**                                                                                             | **File type** | **Description**                                                                                                                       |
|------------------------------------------------------------------------------------------------------|-|---------------------------------------------------------------------------------------------------------------------------------------|
| **Alignment**                                                                                        | :~~: | ~~                                                                                                                                    |
| **`{sample}_{type}_{analysisID}/data/{sample}_{type}_{analysisID}.bam`**                             | bam | Alignment file used for downstream analysis                                                                                           |
| `{sample}_{type}_{analysisID}/data/{sample}_{type}_{analysisID}.bam.bai`                             | bai | Alignment index file                                                                                                                  |
| **SNV and INDELs**                                                                                   | :~~: | ~~                                                                                                                                    |
| `{sample}_{type}_{analysisID}/data/`<br />`{sample}_{type}_{analysisID}.pathology.svdb_query.vcf.gz` | vcf | Final hard-and-soft filtered and annotated file with SNV, INDEL, and CNV? <br /> **OBS! Do not use as input to QCI**                  |
| **`{sample}_{type}_{analysisID}/data/`<br />`{sample}_{type}_{analysisID}.snv-indels.vcf.gz`**       | vcf | Final hard-and-soft filtered and annotated file with SNV and INDEL variants <br />used as input to QCI and to create th Excel report. |


## Reports output files
These files are used by:
- the geneticists to review and analyze the variants that were called,
- the scientists in the wet lab to assess the quality of the sequencing experiment and to validate this.
<br />

| **File**                                                                                         | **File type** | **Description**                                                                                                             |
|--------------------------------------------------------------------------------------------------|---------------|-----------------------------------------------------------------------------------------------------------------------------|
| **SNV and INDELs**                                                                               | :~~:          | ~~                                                                                                                          |
| `{sample}_{type}_{analysisID}/reports/{sample}_{type}_{analysisID}.report.xlsx`                  | xlsx          | Excel report with a compilation of the called variants <br />from different callers. See "Excel" page of the documentation. |
| **CNVs**                                                                                         | :~~:          | ~~                                                                                                                          |
| **`{sample}_{type}_{analysisID}/reports/{sample}_{type}_{analysisID}.pathology.cnv.html`**       | html          | Interactive html report of the filtered CNV results <br />using tumor content estimated by a pathologist                    | |
| **QC**                                                                                           | :~~:          | ~~                                                                                                                          |
| `multiqc_{analysisID}.html`                                                                      | html          | MultiQC report for all samples in the multiplexed sequence run.                                                             |
| **Bamsnap**                                                                                      | :~~:          | ~~                                                                                                                          |
| `{sample}_{type}_{analysisID}/reports/{sample}_{type}_{analysisID}_bamsnap/variant_list.html`    | html          | Report with the IGV screenschots for SNVs with VAF > 5%. <br />To see a screenshot, click on "CHROM" or "POS" in a row      |
| `{sample}_{type}_{analysisID}/reports/{sample}_{type}_{analysisID}_bamsnap/bamsnap_images/*.png` | png           | Standalone screenshots for the SNVs as PNG images.                                                                          |


## Structure of the folder with the output files

The folder with the output files is essentially the home folder `{analysisID}` in which Snakemake is run.
If the pipeline is run without the option `--notemp` in Snakemake, most of the files created during the analysis are deleted upon completion of the pipeline.
Else, the structure of the home folder should be the following one:

```bash
<analysisID>/
├── alignment
│   ├── bwa_mem
│   ├── picard_mark_duplicates
│   ├── samtools_extract_reads
│   └── samtools_merge_bam
├── bamsnap
│   ├── bamsnap
│   │   ├── ...
│   │   └── <sampleID>_<type>
│   │       ├── bamsnap_images
│   │       ├── css
│   │       ├── js
│   │       ├── sample_list
│   │       └── variant_list
│   └── create_pos_list
├── cnv_sv
│   ├── cnvkit_batch
│   │   ├── ...
│   │   └── <sampleID>
│   ├── cnvkit_call
│   ├── cnvkit_vcf
│   ├── gatk_collect_allelic_counts
│   ├── gatk_collect_read_counts
│   ├── gatk_denoise_read_counts
│   ├── gatk_model_segments
│   ├── gatk_vcf
│   ├── pindel
│   ├── pindel_vcf
│   ├── svdb_merge
│   └── svdb_query
├── logs
├── prealignment
│   └── fastp_pe
├── qc
│   ├── fastqc
│   ├── mosdepth_bed_coding
│   ├── multiqc
│   │   └── multiqc_DNA_data
│   ├── picard_collect_alignment_summary_metrics
│   ├── picard_collect_duplication_metrics
│   ├── picard_collect_gc_bias_metrics
│   ├── picard_collect_hs_metrics
│   ├── picard_collect_insert_size_metrics
│   └── samtools_stats
├── reports
│   ├── cnv_html_report
│   └── xlsx
├── results
│   ├── ...
│   └── <sampleID>_<analysisID>
│       ├── data
│       └── reports
│           └── <sampleID>_<type>_<analysisID>_bamsnap
│               ├── bamsnap_images
│               ├── css
│               ├── js
│               ├── sample_list
│               └── variant_list
├── (slurm) # if the pipeline is run on a system that uses the Slurm workload manager
├── snv_indels
│   ├── bcbio_variation_recall_ensemble
│   │   ├── ...
│   │   └── <sampleID>_<type>.ensembled-work
│   ├── bcftools_concat
│   ├── bed_split
│   ├── gatk_mutect2
│   └── vardict
└── versions
    └── software_<run_date>
```

Some extra folders with the source code and configurations for the pipeline might be present as well, 
depending on the setup for the execution.