# Alignment

The steps are executed via the pipeline Poppy GMS. For each sample, they consist in:
* aligning the reads with BWA-mem, 
* splitting the alignments per chromosome for better computational performance in downstream processing,
* marking the duplicate reads, 
* merging the per-chromosome aligned files into a single BAM file.

## Implementation
See the [alignment hydra-genetics module](https://hydra-genetics-alignment.readthedocs.io/en/latest/) 
documentation for more details on the softwares that are used to trim and merge the reads. 
Default hydra-genetics settings/resources are used if no configuration is specified.

## Software versions
See what **containers** are specified in the [config YAML file of Poppy GMS](https://github.com/genomic-medicine-sweden/poppy/blob/v0.2.0/config/config.yaml),
which are possibly overwritten and/or complemented by containers in [poppy-uppsala](https://github.com/clinical-genomics-uppsala/poppy_uppsala/blob/main/config/config_uppsala_nextseq.yaml).

