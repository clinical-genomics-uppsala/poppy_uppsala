# Output files and delivered files

Many temporary files, only a subset of them is forwarded to the geneticists.
These final output files are listed and shortly described in the table below.
`{sample}` is replaced with the sample name and `{type}` is replaced with the sample type (T/N/R, where T=Tumor, N=Normal, R=RNA). 

The **files that are delivered to the geneticists** are shown in bold.

## Structure of the folder with the delivered files

For each sample in the analysis, the following files are delivered:

```bash
<sampleID>_<analysisID>/
  data/
    <sampleID>_<analysisID>.bam
    <sampleID>_<analysisID>.bam.bai
    <sampleID>_<analysisID>.pathology.svdb_query.vcf.gz
    <sampleID>_<analysisID>.snv-indels.vcf.gz
  reports/
multiqc_<analysisID>.html
```

## Reports output files
These files are used by:
- the geneticists to review and analyze the variants that were called,
- the scientists in the wet lab to assess the quality of the sequencing experiment and to validate this.
<br />

| **File**                                                                                              | **File type** | **Description**                                                                                                       |
|-------------------------------------------------------------------------------------------------------|-|-----------------------------------------------------------------------------------------------------------------------|
| **SNV and INDELs**                                                                                    | :~~: | ~~                                                                                                                    |
| `reports/{sample}_{type}_{analysisID}.report.xlsx`                                                    | vcf | Excel report with a compilation of the called variants from different callers. See "Excel" page of the documentation. |
| **CNVs**                                                                                              | :~~: | ~~                                                                                                                    |
| **`reports/{sample}_{type}_{analysisID}.pathology.cnv.html`**                                         | html | Interactive html report of the filtered CNV results <br />using tumor content estimated by a pathologist              |
| `reports/{sample}_{type}_{analysisID}.purecn.cnv.html`                                                | html | Interactive html report of the filtered CNV results <br />using tumor content estimated by PureCN                     |
| **QC**                                                                                                | :~~: | ~~                                                                                                                    |
| `...`                                                                                                 | table | Copied to another location                                                                                            |
| **Bamsnap**                                                                                           | :~~: | ~~                                                                                                                    |
| `reports/{sample}_{type}_{analysisID}_bamsnap`                                                        | table | ...                                                                                                                   |
| `reports/{sample}_{type}_{analysisID}_bamsnap`                                                        | table | ...                                                                                                                   |

## Data output files
These files are used by the geneticists whenever some extra checks are needed under the analysis of the variants.
<br />

| **File**                                                                       | **File type** | **Description** |
|--------------------------------------------------------------------------------|-|-|
| **Alignment**                                                                  | :~~: | ~~ |
| **`data/{sample}_{type}_{analysisID}.bam`**                                    | bam | Alignment file used for downstream analysis|
| `data/{sample}_{type}_{analysisID}.bam.bai`                                    | bai | Alignment index file |
| **SNV and INDELs**                                                             | :~~: | ~~ |
| **`data/`<br />`{sample}_{type}_{analysisID}.pathology.svdb_query.vcf.gz`**    | vcf | Final hard filtered and annotated file with SNV and INDEL variants. <br /> **OBS! Do not use as input to QCI** |
| **`data/`<br />`{sample}_{type}_{analysisID}.snv-indels.vcf.gz`**              | vcf | Final hard filtered and annotated file with SNV and INDEL variants <br />used as input to QCI |

