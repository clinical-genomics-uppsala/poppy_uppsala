# Softwares used in the poppy_uppsala module

[//]: # ()
[//]: # (## [results_report_bedtools_intersect]&#40;https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html&#41;)

[//]: # (Bedtools intersect to extract lowcoverage regions inside coding exons.)

[//]: # ()
[//]: # (### :snake: Rule)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_SOURCE__results_report__results_report_bedtools_intersect#)

[//]: # ()
[//]: # (#### :left_right_arrow: input / output files)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_TABLE__results_report__results_report_bedtools_intersect#)

[//]: # ()
[//]: # (### :wrench: Configuration)

[//]: # ()
[//]: # (#### Software settings &#40;`config.yaml`&#41;)

[//]: # ()
[//]: # (#CONFIGSCHEMA__results_report_bedtools_intersect#)

[//]: # ()
[//]: # (#### Resources settings &#40;`resources.yaml`&#41;)

[//]: # ()
[//]: # (#RESOURCESSCHEMA__results_report_bedtools_intersect#)

[//]: # ()
[//]: # ()
[//]: # (## results_report_xlsx)

[//]: # (python script that summerizes results into xlsx report)

[//]: # ()
[//]: # (### :snake: Rule)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_SOURCE__results_report__results_report_xlsx#)

[//]: # ()
[//]: # (#### :left_right_arrow: input / output files)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_TABLE__results_report__results_report_xlsx#)

[//]: # ()
[//]: # (### :wrench: Configuration)

[//]: # ()
[//]: # (#### Software settings &#40;`config.yaml`&#41;)

[//]: # ()
[//]: # (#CONFIGSCHEMA__results_report_xlsx#)

[//]: # ()
[//]: # (#### Resources settings &#40;`resources.yaml`&#41;)

[//]: # ()
[//]: # (#RESOURCESSCHEMA__results_report_xlsx#)

[//]: # ()
[//]: # ()
[//]: # (## bamsnap_create_pos_list)

[//]: # (A python script to generate bedfile with PASS and AF >= params.af for bamsnap-rule)

[//]: # ()
[//]: # (### :snake: Rule)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_SOURCE__bamsnap__bamsnap_create_pos_list#)

[//]: # ()
[//]: # (#### :left_right_arrow: input / output files)

[//]: # ()
[//]: # (#SNAKEMAKE_RULE_TABLE__bamsnap__bamsnap_create_pos_list#)

[//]: # ()
[//]: # (### :wrench: Configuration)

[//]: # ()
[//]: # (#### Software settings &#40;`config.yaml`&#41;)

[//]: # ()
[//]: # (#CONFIGSCHEMA__bamsnap_create_pos_list#)

[//]: # ()
[//]: # (#### Resources settings &#40;`resources.yaml`&#41;)

[//]: # ()
[//]: # (#RESOURCESSCHEMA__bamsnap_create_pos_list#)


## [bamsnap](https://bamsnap.readthedocs.io/en/latest/)
To automatically generate snapshots of variants bamsnap is used

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bamsnap__bamsnap#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bamsnap__bamsnap#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bamsnap#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bamsnap#


## bamsnap_hd829
Dummy rule for HD829 samples instead of bamsnap rule

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bamsnap__bamsnap_hd829#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bamsnap__bamsnap_hd829#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bamsnap_hd829#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bamsnap_hd829#

## version_update_poppy
Since pipeline version can only use current git repo, the poppy-version-file.yaml will contain poppy_uppsala version, this simple shell script overwrite the file with the actual poppy-version taken from config.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__version__version_update_poppy#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__version__version_update_poppy#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__version_update_poppy#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__version_update_poppy#
