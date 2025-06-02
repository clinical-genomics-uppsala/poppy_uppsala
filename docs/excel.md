# Description of the sheets in the Excel report

One Excel report is created for each sample in the analysis.

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

All filters can be reset and customized.

Usual coverage range over the exons: ca. 700-800x if NextSeq machine was used for sequencing, ca. 1200-1400x with NovaSeqX.

Some regions are systematically reported in the "Low Coverage" tab. 
For instance the exon 1 in GNAS (GNAS[1]) is systematically listed if the sequencing was done on NextSeq.
However, using NovaSeqX ensures sufficient coverage in that region.

## SNVs sreenshots from bamsnap

The images are not pasted in the report but they are available in another directory.