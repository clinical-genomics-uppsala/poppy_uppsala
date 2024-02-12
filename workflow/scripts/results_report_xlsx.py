#!/bin/python3

from results_report_create_tables import *
import sys
from datetime import date

vcf = snakemake.input.vcf
pindel = snakemake.input.pindel
sequenceid = snakemake.params.sequenceid

sample = snakemake.params.sample
""" Create data tables to populate excel sheets with """
snv_table = create_snv_table(vcf, sequenceid)
pindel_table = create_pindel_table(pindel, sequenceid)
known_table, known_percent = create_known_variants_table(vcf, pindel, sequenceid)

pindel_genes = ""  # To be done

""" xlsx file with sheets """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
if sample.lower() == "hd829":
    worksheet_known = workbook.add_worksheet("Known variants")
worksheet_snv = workbook.add_worksheet("SNVs")
worksheet_pindel = workbook.add_worksheet("Pindel")


empty_list = ["", "", "", "", "", ""]
format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_line = workbook.add_format({"top": 1})
format_orange = workbook.add_format({"bg_color": "#ffd280"})
format_red = workbook.add_format({"font_color": "red"})
format_table_heading = workbook.add_format({"bold": True, "text_wrap": True})

# Overview
worksheet_overview.write(0, 0, sample, format_heading)
worksheet_overview.write(1, 0, "RunID: " + sequenceid)
worksheet_overview.write(2, 0, "Processing date: " + date.today().strftime("%B %d, %Y"))
worksheet_overview.write_row(3, 0, empty_list, format_line)

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")
worksheet_overview.write_row(6, 0, empty_list, format_line)

worksheet_overview.write(7, 0, "Sheets:", format_table_heading)
i = 8
if sample.lower() == "hd829":
    worksheet_overview.write_url(i, 0, "internal:'Known variants'!A1", string="Known variants")
    i += 1
worksheet_overview.write_url(i, 0, "internal:'SNVs'!A1", string="SNVs identified")
worksheet_overview.write_url(i + 1, 0, "internal:'Pindel'!A1", string="Pindel results")

if sample.lower() == "hd829":
    worksheet_overview.write(i + 3, 0, "Percent of known variants identified:")
    if known_percent < 1:
        worksheet_overview.write(i + 4, 0, str(known_percent * 100) + " %", format_red)
    else:
        worksheet_overview.write(i + 4, 0, str(known_percent * 100) + " %")
    i += 4

i = i + 2
worksheet_overview.write(i, 0, "Artifact panel used: " + snakemake.params.artifact)
# bedfiles?
worksheet_overview.write(i + 1, 0, "Twist Myeloid capture panel TE-91316667_hg38 used")
worksheet_overview.write(
    i + 2,
    0,
    "with the pipeline Poppy (v"
    + snakemake.params.poppy_version
    + ") and local uppsala implementations (v"
    + snakemake.params.uppsala_version
    + ").",
)
worksheet_overview.write_url(i + 4, 0, "https://gms-poppy.readthedocs.io/en/latest/", string="Poppy documentation")
worksheet_overview.write_url(
    i + 5, 0, "https://github.com/clinical-genomics-uppsala/poppy_uppsala/tree/develop", string="Uppsala configurations"
)

# SNVs sheet
for x in VariantFile(vcf).header.records:
    if x.key == "VEP":
        vep_line = x.value

worksheet_snv.set_column(5, 5, 10)
worksheet_snv.set_column(11, 13, 10)
worksheet_snv.write("A1", "Variants found", format_heading)
worksheet_snv.write("A3", "Sample: " + str(sample))
worksheet_snv.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_snv.write("A6", "Databases used: " + vep_line)

worksheet_snv.write("A8", "Filters: ", format_orange)
worksheet_snv.write("B9", "DP_200: Soft filter on depth lower than 200X", format_orange)
worksheet_snv.write("B10", "AD_5: Soft filter variants with few observations (AD lower than 5)", format_orange)
worksheet_snv.write("B11", "Artifact_gt_3: Soft filter variants found in more than 3 normal samples", format_orange)
worksheet_snv.write(
    "B12", "PopAF_0.02: : Soft filter germline if >2 % in any population from 1000 genomes, ESP or gnomADe", format_orange
)
worksheet_snv.write(
    "B13",
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic" +
    " ID on the position",
    format_orange,
)

table_area = "A18:V" + str(len(snv_table["data"]) + 18)

worksheet_snv.add_table(table_area, {"data": snv_table["data"], "columns": snv_table["headers"], "style": "Table Style Light 1"})
table_area_data = "A19:V" + str(len(snv_table["data"]) + 19)
cond_formula = '=LEFT($A19, 4)<>"PASS"'
worksheet_snv.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

# Pindel sheet
worksheet_pindel.set_column(5, 5, 10)
worksheet_pindel.set_column(11, 13, 10)
worksheet_pindel.write("A1", "Variants found", format_heading)
worksheet_pindel.write("A3", "Sample: " + str(sample))
worksheet_pindel.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_pindel.write("A6", "To limit runtime pindel were used with the designfile used: " + snakemake.params.pindel_bed)
worksheet_pindel.write("B7", "Which includes the following genes: ")
i = 7
for gene in pindel_genes:
    worksheet_pindel.write("C" + str(i), gene)
    i += 1

worksheet_pindel.write("A" + str(i + 1), "Filters: ", format_orange)
worksheet_pindel.write("B" + str(i + 2), "DP_200: Soft filter on depth lower than 200X", format_orange)
worksheet_pindel.write("B" + str(i + 3), "AD_5: Soft filter variants with few observations (AD lower than 5)", format_orange)
worksheet_pindel.write("B" + str(i + 4), "Artifact_gt_3: Soft filter variants found in more than 3 normal samples", format_orange)
worksheet_pindel.write(
    "B" + str(i + 5),
    "PopAF_0.02: Soft filter germline if >2 % in any population from 1000 genomes, ESP or gnomADe",
    format_orange,
)
worksheet_pindel.write(
    "B" + str(i + 6),
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic " +
    "ID on the position",
    format_orange,
)
i = i + 8
table_area = "A" + str(i) + ":T" + str(len(pindel_table["data"]) + i)
worksheet_pindel.add_table(
    table_area, {"data": pindel_table["data"], "columns": pindel_table["headers"], "style": "Table Style Light 1"}
)

table_area_data = "A" + str(i + 1) + ":T" + str(len(pindel_table["data"]) + i + 1)
cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
worksheet_pindel.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

if sample.lower() == "hd829":
    worksheet_known.set_column(3, 4, 10)
    worksheet_known.set_column(6, 6, 10)
    worksheet_known.write("A1", "Variants known for HD829", format_heading)
    worksheet_known.write("A3", "Sample: " + str(sample))
    i = 5
    table_area = "A" + str(i) + ":O" + str(len(known_table["data"]) + i)
    worksheet_known.add_table(
        table_area, {"data": known_table["data"], "columns": known_table["headers"], "style": "Table Style Light 1"}
    )

# # coverage
workbook.close()
