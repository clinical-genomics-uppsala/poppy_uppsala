# Output files and delivered files

Many temporary files, only a subset of them is forwarded to the geneticists.
These final output files are listed and shortly described in the table below.
`{sample}` is replaced with the sample name and `{type}` is replaced with the sample type (T/N/R, where T=Tumor, N=Normal, R=RNA). 
The most important files and shown in bold.

## Structure of the folder with the delivered files

```bash
results/
  data/
  reports/
```

## Reports output files
These files are used by:
- the geneticists to review and analyze the variants that were called,
- the scientists in the wet lab to assess the quality of the sequencing experiment and to validate this.
<br />

| **File**                                                                                              | **File type** | **Description**                                                                                 |
|-------------------------------------------------------------------------------------------------------|-|-------------------------------------------------------------------------------------------------|
| **SNV and INDELs**                                                                                    | :~~: | ~~                                                                                              |
| `results/dna/additional_files/vcf/{caller}_{sample}_{type}.vcf.gz`                                    | vcf | SNV and INDEL variants file from the individual callers                                         |
| `results/dna/additional_files/vcf/{sample}_{type}.annotated.vcf.gz`                                   | vcf | Merged and annotated SNV and INDEL variants file                                                |
| `results/dna/additional_files/vcf/`<br />`{sample}_{type}.annotated.exon_only.filter.soft_filter.vcf` | vcf | Soft filtered, merged and annotated SNV and INDEL variants file <br />without double variants   |
| `results/dna/additional_files/vcf/`<br />`{sample}_{type}.annotated.exon_only.filter.hard_filter.vcf` | vcf | Hard filtered, merged and annotated SNV and INDEL variants file <br />without double variants   |
| **CNVs**                                                                                              | :~~: | ~~                                                                                              |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.manta_tumorSV.vcf.gz`               | vcf | Manta variant calling                                                                           |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.cnvkit.scatter.png`                 | image | CNVkit genome CNV plot                                                                          |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.cnvkit.diagram.pdf`                 | image | CNVkit chromosome CNV plot                                                                      |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.purecn.svdb_query.vcf`              | vcf | Merged CNV vcf by SVDB from the two callers <br />using tumor content estimated by PureCN       |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.pathology.svdb_query.vcf`           | vcf | Merged CNV vcf by SVDB from the two callers <br />using tumor content estimated by a pathologist |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.deletions.tsv`                      | table | Table of called small CNVs deletions by in-house script                                         |
| `results/dna/additional_files/cnv/{sample}_{type}/{sample}_{type}.amplifications.tsv`                 | table | Table of called small CNVs amplifications by in-house script                                    |
| **SVs**                                                                                               | :~~: | ~~                                                                                              |
| `results/dna/additional_files/fusion/{sample}_{type}.gene_fuse_fusions.txt`                           | table | ...                                                                                             |
| `results/rna/additional_files/fusion/{sample}_{type}.arriba.fusions.tsv`                              | table | ...                                                                                             |
| `results/rna/additional_files/fusion/{sample}_{type}.star-fusion.fusion_predictions.tsv`              | table | ...                                                                                             |
| `results/rna/additional_files/fusion/{sample}_{type}.fusioncatcher.fusion_predictions.txt`            | table | ...                                                                                             |
| **QC**                                                                                                | :~~: | ~~                                                                                              |
| `results/dna/additional_files/qc/{sample}_{type}.alignment_summary_metrics.txt`                       | table | ...                                                                                             |
| `results/dna/additional_files/qc/{sample}_{type}.contamination.table`                                 | table | ...                                                                                             |
| `results/dna/additional_files/qc/{sample}_{type}.duplication_metrics.txt`                             | table | ...                                                                                             |
| **Bamsnap**                                                                                           | :~~: | ~~                                                                                              |
| `results/dna/additional_files/qc/{sample}_{type}.insert_size_metrics.txt`                             | table | ...                                                                                             |
| `results/dna/additional_files/qc/{sample}_{type}.samtools-stats.txt`                                  | table | ...                                                                                             |

## Data output files
These files are used by the geneticists whenever some extra checks are needed under the analysis of the variants.
<br />

| **File**                                                                                               | **File type** | **Description** |
|--------------------------------------------------------------------------------------------------------|-|-|
| **Alignment**                                                                                          | :~~: | ~~ |
| **`bam_dna/{sample}_{type}.bam`**                                                                      | bam | Alignment file used for downstream analysis|
| `bam_dna/{sample}_{type}.bam.bai`                                                                      | bai | Alignment index file |
| **`bam_rna/{sample}_{type}.star_fusion.bam`**                                                          | bam | Alignment file  |
| `bam_rna/{sample}_{type}.star_fusion.bam.bai`                                                          | bai | Alignment index file |
| `bam_dna/mutect2_indel_bam/{sample}_{type}.bam`                                                        | bam | Realigned regions around INDELs, used to look at indels in for example IGV |
| `bam_dna/mutect2_indel_bam/{sample}_{type}.bam.bai`                                                    | bai | Alignment index file |
| **SNV and INDELs**                                                                                     | :~~: | ~~ |
| **`results/dna/vcf/`<br />`{sample}_{type}.annotated.exon_only.filter.hard_filter.codon_snv.vcf`**     | vcf | Final hard filtered and annotated file with SNV and INDEL variants. <br /> **OBS! Do not use as input to QCI** |
| **`results/dna/vcf/`<br />`{sample}_{type}.annotated.exon_only.filter.hard_filter.codon_snv.qci.vcf`** | vcf | Final hard filtered and annotated file with SNV and INDEL variants <br />used as input to QCI |
| `gvcf_dna/{sample}_{type}.mosdepth.g.vcf.gz`                                                           | g.vcf | Genomic vcf file with coverage in all position in the design |
| **CNVs**                                                                                               | :~~: | ~~ |
| **`results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv.html`**                               | html | Interactive html report of the filtered CNV results <br />using tumor content estimated by a pathologist |
| `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv.html`                                      | html | Interactive html report of the filtered CNV results <br />using tumor content estimated by PureCN |
| **`results/dna/cnv/{sample}_{type}/{sample}_{type}.pathology.cnv_report.tsv`**                         | table | Table of filtered CNVs using tumor content <br />estimated by a pathologist |
| `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.cnv_report.tsv`                                | table | Table of filtered CNVs using tumor content <br />estimated by PureCN |
| `results/dna/cnv/{sample}_{type}/{sample}_{type}.purecn.purity.txt`                                    | table | Estimated tumor content and ploidity by PureCN |
| **SVs**                                                                                                | :~~: | ~~ |
| **`results/dna/fusion/{sample}_{type}.gene_fuse_report.tsv`**                                          | table | Filtered report of DNA fusions from GeneFuse |
| **`results/rna/fusion/{sample}_{type}.fusion_report.tsv`**                                             | table | Filtered report of RNA fusions from three different callers |
| **`results/rna/fusion/{sample}_{type}.exon_skipping.tsv`**                                             | table | Called exon skipping in MET and EGFR |
| `results/rna/fusion/{sample}_{type}.arriba.fusions.pdf`                                                | image | Image generated by Arriba visualizing the called fusions |

