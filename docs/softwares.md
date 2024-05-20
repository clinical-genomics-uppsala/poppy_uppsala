# Softwares used in the poppy_uppsala module

## [results_report_bedtools_intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
Bedtools intersect to extract lowcoverage regions inside coding exons.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__results_report__results_report_bedtools_intersect#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__results_report__results_report_bedtools_intersect#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__results_report_bedtools_intersect#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__results_report_bedtools_intersect#


## results_report_xlsx
python script that summerizes results into xlsx report

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__results_report__results_report_xlsx#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__results_report__results_report_xlsx#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__results_report_xlsx#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__results_report_xlsx#


## bamsnap_create_pos_list
A python script to generate bedfile with PASS and AF >= params.af for bamsnap-rule

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bamsnap__bamsnap_create_pos_list#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bamsnap__bamsnap_create_pos_list#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bamsnap_create_pos_list#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bamsnap_create_pos_list#


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