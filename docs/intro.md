# Introduction

An introduction to poppy-uppsala.
Bioinformatics for the analysis panel of genes Twist Myeloid.

## The Twist Myeloid gene panel
The **Twist Myeloid gene panel** is a gDNA target capture panel for comprehensive genomic profiling of leukemias.
The panel was designed and validated as a collaborative effort between hematology teams within Genomic Medicine Sweden (GMS). 
Reference genome GRCh38.
Replaces pomfrey GRCh37.
The design allows for clinical use as well as for research applications.

Typically, the sequencing experiments are performed in multiplex settings on machines like NextSeq and NovaSeqX.

## Pipeline
The poppy pipeline analyses tumor-only blood/bone marrow samples and reports small variants, CNVs against the Twist Myeloid gene panel.
The pipeline is build using the hydra-genetics framework. 
All new releases are run through real data sets and compared to previous versions to ensure that the expected result are generated.