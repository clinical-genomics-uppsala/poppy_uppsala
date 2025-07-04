# Excel report

One Excel report is created for each sample in the analysis and generated from the VCF files that are produced by the 
pipeline.
The report is named `<sampleID>_<type>_<analysisID>_report.xlsx` 
and contains several sheets with information about the variants found in the sample, 
as well as coverage metrics.
The variants reported in the Excel file are filtered and curated to facilitate the interpretation
of the genetic profile of the patients against the panel of genes "Twist Myeloid".

The script that generates the Excel report is `poppy_uppsala/scripts/excel_report.py` and it is based on the Python
package `XlsxWriter`.

## General filtering mechanism

The Excel report contains several sheets with different types of variants.
For each variant type, the filtering strategy is based on:
- hard filters: they are applied to the VCF files before the Excel report is generated,
- soft filters: additional filters that are applied in the Excel report itself.

The filtering criteria can be customized by the user via the configuration files in the `config/filters` 
directory. 
Overall, we have used mostly the VEP annotations and the VAF to filter the variants.

## Description of the sheets in the Excel report

| **Tab**      | **Description**                                                                                                                                                                                 |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Overview     | Information about the BED files that were used in the analysis. Hyperlinks to the other sheets in the report and to the screenshots generated with bamsnap. Summary of the coverage thresholds. |
| CLL          | List of variants that are relevant to the diagnostic of CLL (subset from "SNVs"). The list is prefiltered and shows only variants with > 2 % AF and filter-flag PASS.                           |
| MYELOID      | List of variants that are relevant to the diagnostic of myeloid leukemias (subset from "SNVs"). The list is prefiltered and shows only variants with > 2 % AF and filter-flag PASS.             |
| SNVs         | List of SNVs and indels found by one or more callers (GATK Mutect2, VarDict. The list is prefiltered and shows only variants with > 2 % AF and filter-flag PASS.                                |
| Pindel       | List of indels found by Pindel. Variant calling with Pindel is performed in a restricted set of regions in order to cut the computational cost of the analysis.                                 |
| Intron       | List of variants identified in a restricted set of intron regions. Some variants might be deemed to another consequence than "intron_variant" though.                                           |
| Synonymous   | List of positions where synonymous variants were found.                                                                                                                                         |
| Low Coverage | Regions with coverage below 100x and the corresponding transcripts.                                                                                                                             |
| Coverage     | Average coverage of each region in exon bedfile with the gene, the exon number, and the relevant transcript.                                                                                    |
| QCI          | ?                                                                                                                                                                                               |

All variants should be annotated with VEP.

Upon opening the Excel report, a default filtering is applied to the variants that are shown in the sheets.
All filters can be reset and customized by clicking on the header of a column in a table and (un)checking boxes to 
show only the variants matching the filtering criteria, for instance all variants on chr2 with VAF > 10 %.
When the default filtering is active, all boxes are checked in all columns.

Usual coverage range over the exons: ca. 700-800x if NextSeq machine was used for sequencing, ca. 1200-1400x with NovaSeqX.

Some regions are systematically reported in the "Low Coverage" tab. 
For instance the exon 1 in GNAS (GNAS[1]) is systematically listed if the sequencing was done on NextSeq.
However, using NovaSeqX ensures sufficient coverage in that region.

## SNVs screenshots from bamsnap

The images are not pasted in the report but they are available in another directory.

## How to customize the report

Modify the filtering criteria in the configuration files in the `config/filters` directory.

We refer the user to the [XlsxWriter documentation](https://xlsxwriter.readthedocs.io/) 
for more details on how to modify the script that generates the report.