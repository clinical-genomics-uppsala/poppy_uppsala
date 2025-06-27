# Prealignment

The steps are executed via the pipeline Poppy GMS, they consist in trimming the reads and 
merging together the fastq files that belong to the same sample.
Typically, if the fastq files derive from a NextSeq machine that ran PE sequencing,
there are 8 fastq files for each sample in the analysis:

* \<analysisID>_\<sampleID>_L001_R1_001.fastq.gz: Read 1 from lane 1	
* \<analysisID>_\<sampleID>_L002_R1_001.fastq.gz: Read 2 from lane 2
* \<analysisID>_\<sampleID>_L003_R1_001.fastq.gz: Read 1 from lane 3	
* \<analysisID>_\<sampleID>_L004_R1_001.fastq.gz: Read 1 from lane 4
* \<analysisID>_\<sampleID>_L001_R2_001.fastq.gz: Read 2 from lane 1	
* \<analysisID>_\<sampleID>_L002_R2_001.fastq.gz: Read 2 from lane 2
* \<analysisID>_\<sampleID>_L003_R2_001.fastq.gz: Read 2 from lane 3	
* \<analysisID>_\<sampleID>_L004_R2_001.fastq.gz: Read 2 from lane 4

## Implementation
See the [prealignment hydra-genetics module](https://hydra-genetics-prealignment.readthedocs.io/en/latest/) 
documentation for more details on the softwares that are used to trim and merge the reads. 
Default hydra-genetics settings/resources are used if no configuration is specified.

## Software versions
See what **containers** are specified in the [config YAML file of Poppy GMS](https://github.com/genomic-medicine-sweden/poppy/blob/v0.2.0/config/config.yaml),
which are possibly overwritten and/or complemented by containers in [poppy-uppsala](https://github.com/clinical-genomics-uppsala/poppy_uppsala/blob/main/config/config_uppsala_nextseq.yaml).

[//]: # (<br />)

[//]: # (![dag plot]&#40;images/prealignment.png&#41;{: style="height:18%;width:18%"})

[//]: # ()
[//]: # (## Pipeline output files:)

[//]: # (Only temporary intermediate files are created.)

[//]: # ()
[//]: # (## Trimming)

[//]: # (Trimming of fastq files is performed by **[fastp]&#40;https://github.com/OpenGene/fastp&#41;** v0.20.1.  )

[//]: # ()
[//]: # (### Configuration)

[//]: # ()
[//]: # ()
[//]: # (**Resources**)

[//]: # ()
[//]: # (| **Options** | **Value** |)

[//]: # (|-------------|-|)

[//]: # (| mem_mb | 30720 |)

[//]: # (| mem_per_cpu | 6144 |)

[//]: # (| threads | 5 |)

[//]: # ()
[//]: # (## Merging)

[//]: # (Merging of fastq files belonging to the same sample are performed by simply concatenating the files with **cat**.)

[//]: # ()
[//]: # (<br />)
