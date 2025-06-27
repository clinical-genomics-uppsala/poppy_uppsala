# Overview of the pipeline
Here is a brief overview of the entire pipeline. For details see subsections and the [hydra-genetics](https://github.com/hydra-genetics/hydra-genetics) documentation.

Essentially, a few extra steps specific to the needs of the geneticists at Clinical Genomics Uppsala are added
on the top of the pipeline [Poppy (Genomic Medicine Sweden)](https://github.com/genomic-medicine-sweden/poppy).


## Main processing steps from poppy GMS

Depending on the sequencing setup that is used, a **preliminary demultiplexing** of the samples may be required.
We assume in the following steps that the fastq files are demultiplexed.

1. **Input files**: directory with fastq files, also used to create the `samples.tsv` and `units.tsv` files (see section "Running Poppy Uppsala")
2. **Trimming** using fastp and **merging**
3. **Alignment** using BWA-mem
4. **Mark duplicates** using Picard
5. **SNV and INDEL**  
  5.1 Calling using Mutect2 and Vardict  
  5.2 Annotation using VEP and hydra-genetics  
  5.3 Filtering using bcftools and hydra-genetics  
6. **CNV**  
  6.1 Calling using CNVkit and GATK CNV  
  6.2 Merging using SVDB  
  6.3 Annotation using SVDB and hydra-genetics  
  6.4 Filtering using hydra-genetics  
  6.5 CNV html report using hydra-genetics
7. **Pindel** for more complex indels  
  8.1 Calling using Pindel in chosen restricted regions  
  8.2 Annotation using VEP and hydra-genetics  
  8.3 Filtering using bcftools and hydra-genetics  
9. **QC**  
  9.1 QC measures from Samtools, Picard, FastQC, GATK  
  9.2 MultiQC hmtl report  
  9.3 Coverage report from Mosdepth  

## Additional steps specific to poppy-uppsala
10. Calculate coverage only in chosen exon regions,
11. Produce an Excel report with filtered and curated variants from different callers (Mutect2, VarDict, Pisces, GATK, Pindel), 
as well as coverage metrics. The report has been designed to address requests from the geneticists' team 
in order to facilitate the interpretation of the genetic profile of the patients against the panel of genes "Twist Myeloid".
12. Produce a graphical report with MultiQC and reorder the samples such that they are displayed in the same order 
as in the samplesheet used in the wet lab.
13. Take an automated screenshot in IGV of the filtered variants that are located in genes in the panel and that have a VAF > 5%.

## Rule graph
![rulegraph](https://raw.githubusercontent.com/clinical-genomics-uppsala/poppy_uppsala/patch-readthedocs/images/rulegraph.png)
