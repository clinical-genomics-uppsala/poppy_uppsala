#!/bin/python3

from results_report_create_tables import *
import sys
from datetime import date
from operator import itemgetter

vcf = snakemake.input.vcf
pindel = snakemake.input.pindel
sequenceid = snakemake.params.sequenceid

sample = snakemake.params.sample
min_cov = 200
wanted_transcripts = []

""" Create data tables to populate excel sheets with """
snv_table = create_snv_table(vcf, sequenceid)
pindel_table = create_pindel_table(pindel, sequenceid)
known_table, known_percent = create_known_variants_table(vcf, pindel, sequenceid)

# List genes in pindel bed
pindel_genes = []
with open(snakemake.params.pindelbed, "r") as pindel_file:
    for lline in pindel_file:
        line = lline.strip().split("\t")
        pindel_genes.append(line[3])
pindel_genes = list(dict.fromkeys(pindel_genes))

# Avg per exon/coding region
regionscov_table = {"data": [], "headers": []}
regionscov_table["headers"] = [
    {"header": "Chr"},
    {"header": "Start"},
    {"header": "Stop"},
    {"header": "Gene"},
    {"header": "Exon"},
    {"header": "Transcript"},
    {"header": "Avg Coverage"},
]
bed_table = []
with gzip.open(snakemake.input.mosdepth_regions, "rt") as regions_file:
    for lline in regions_file:
        line = lline.strip().split("\t")
        gene = line[3].split("_")[0]
        transcript = "_".join(line[3].split("_")[1:3])
        exon = str(line[3].split("_")[3])
        coverage_row = [line[0], line[1], line[2], gene, exon, transcript, float(line[4])]
        if coverage_row not in regionscov_table["data"]:
            regionscov_table["data"].append(coverage_row)
        if line[0:5] not in bed_table:
            bed_table.append(line[0:5])

# Low coverage data
lowcov_lines = []
with open(snakemake.input.mosdepth_perbase, "r") as mosdepth_perbase:
    for lline in mosdepth_perbase:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(min_cov) and line[0:4] not in lowcov_lines:
            lowcov_lines.append(line[0:4])
lowcov_lines = sorted(lowcov_lines, key=itemgetter(0, 1))  # Sort based on chr and start pos

lowcov_table = {"data": [], "headers": []}
lowcov_table["headers"] = [
    {"header": "Chr"},
    {"header": "Start"},
    {"header": "Stop"},
    {"header": "Mean Coverage"},
    {"header": "Preferred transcript"},
    {"header": "All transcripts"},
]

num_low_regions = 0
for line in lowcov_lines:
    line[3] = int(line[3])
    exons = []
    for bed_line in bed_table:
        # get all exons that cover that low cov line
        if line[0] == bed_line[0] and int(line[1]) >= int(bed_line[1]) and int(line[2]) <= int(bed_line[2]):
            exons.append(bed_line[3])

    if len(exons) > 0:
    #     if any(exon in wanted_transcripts for exon in exons):
    #         lowcov_table["data"].append(line + list(set(exons) & set(wanted_transcripts)) + [";".join(exons)])
    #         num_low_regions += 1
    #     else:
    #         lowcov_table["data"].append(line + [""] + [";".join(exons)])
        lowcov_table["data"].append(line + [""] + [";".join(exons)])



""" xlsx file with sheets """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
if sample.lower() == "hd829":
    worksheet_known = workbook.add_worksheet("Known variants")
worksheet_snv = workbook.add_worksheet("SNVs")
worksheet_pindel = workbook.add_worksheet("Pindel")
worksheet_lowcov = workbook.add_worksheet("Low Coverage")
worksheet_cov = workbook.add_worksheet("Coverage")

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
worksheet_overview.write_url(i + 2, 0, "internal:'Low Coverage'!A1", string="Low Coverage regions")
worksheet_overview.write_url(i + 3, 0, "internal:'Coverage'!A1", string="Coverage")

if sample.lower() == "hd829":
    i += 5
    worksheet_overview.write(i, 0, "Percent of known variants identified:")
    if known_percent < 1:
        worksheet_overview.write(i + 1, 0, str(known_percent * 100) + " %", format_red)
    else:
        worksheet_overview.write(i + 1, 0, str(known_percent * 100) + " %")
i += 5
worksheet_overview.write(i, 0, "Number of regions in wanted transcripts not coverage by at least " + str(min_cov) + "x: ")
worksheet_overview.write(i+1, 0, str(num_low_regions))

i += 4
worksheet_overview.write(i, 0, "Full design bedfile: " + snakemake.params.bedfile)
worksheet_overview.write(i + 1, 0, "Coding exons bedfile: " + snakemake.params.exonbed)
worksheet_overview.write(i + 2, 0, "Artifact panel used: " + snakemake.params.artifact)

i += 4
worksheet_overview.write(i, 0, "Twist Myeloid capture panel TE-91316667_hg38 used")
worksheet_overview.write(
    i + 1,
    0,
    "with the pipeline Poppy (v"
    + snakemake.params.poppy_version
    + ") and local uppsala implementations (v"
    + snakemake.params.uppsala_version
    + ").",
)
worksheet_overview.write_url(i + 3, 0, "https://gms-poppy.readthedocs.io/en/latest/", string="Poppy documentation")
worksheet_overview.write_url(
    i + 4, 0, "https://github.com/clinical-genomics-uppsala/poppy_uppsala/tree/develop", string="Uppsala configurations"
)

""" SNVs sheet """
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
    "B12", "PopAF_0.02: Soft filter germline if >2 % in any population from 1000 genomes, ESP or gnomADe", format_orange
)
worksheet_snv.write(
    "B13",
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic"
    + " ID on the position",
    format_orange,
)
worksheet_snv.write(
    "B14",
    "Consequence: Soft filter variants which consequence is deemed irrelevant " +
    "(intergenic_variant, NMD_transcript_variant, non_coding_transcript_variant, upstream_gene_variant," +
    " downstream_gene_variant, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, " +
    "regulatory_region_amplification, regulatory_region_variant)",
    format_orange,
)
worksheet_snv.write("B15", "Biotype: Soft filter variants not annotated as protein_coding", format_orange)

worksheet_snv.write("A17", "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. You can then use the drop-downs in the header row to filter to your liking.")
table_area = "A19:V" + str(len(snv_table["data"]) + 19)

worksheet_snv.add_table(table_area, {"columns": snv_table["headers"], "style": "Table Style Light 1"})
table_area_data = "A20:V" + str(len(snv_table["data"]) + 20)
cond_formula = '=LEFT($A20, 4)<>"PASS"'
worksheet_snv.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

# Hide lines that are not PASS by default
i=19
worksheet_snv.autofilter(table_area)
worksheet_snv.filter_column('A', 'Filter != PASS')
for row_data in snv_table["data"]:
    if row_data[0] == "PASS":
        pass
    else:
        worksheet_snv.set_row(i, options={'hidden': True})
    worksheet_snv.write_row(i, 0, row_data)
    i += 1


""" Pindel sheet """
worksheet_pindel.set_column(5, 5, 10)
worksheet_pindel.set_column(11, 13, 10)
worksheet_pindel.write("A1", "Variants found", format_heading)
worksheet_pindel.write("A3", "Sample: " + str(sample))
worksheet_pindel.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_pindel.write("A6", "To limit runtime pindel were used with the designfile used: " + snakemake.params.pindelbed)
worksheet_pindel.write("B7", "Which includes the following genes: ")
i = 8
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
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic "
    + "ID on the position",
    format_orange,
)
worksheet_pindel.write("A"+str(i+8), "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. You can then use the drop-downs in the header row to filter to your liking.")

i += 10
table_area = "A" + str(i) + ":T" + str(len(pindel_table["data"]) + i + 1)
worksheet_pindel.add_table(table_area, {"columns": pindel_table["headers"], "style": "Table Style Light 1"})

table_area_data = "A" + str(i + 1) + ":T" + str(len(pindel_table["data"]) + i)
cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
worksheet_pindel.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

i+=1
worksheet_pindel.autofilter(table_area)
worksheet_pindel.filter_column('A', 'Filter != PASS')
for row_data in pindel_table["data"]:
    if row_data[0] == "PASS":
        pass
    else:
        worksheet_pindel.set_row(i, options={'hidden': True})
    worksheet_pindel.write_row(i, 0, row_data)
    i += 1

""" Known sheet """
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

""" Low coverage """
worksheet_lowcov.set_column(1, 2, 10)
worksheet_lowcov.set_column(4, 5, 25)

worksheet_lowcov.write(0, 0, "Mosdepth low coverage analysis", format_heading)
worksheet_lowcov.write_row(1, 0, empty_list, format_line)
worksheet_lowcov.write(2, 0, "Sample: " + str(sample))
worksheet_lowcov.write(3, 0, "Gene regions with coverage lower than " + str(min_cov) + "x.")

table_area = "A6:F" + str(len(lowcov_table["data"]) + 6)
worksheet_lowcov.add_table(table_area, {"data": lowcov_table["data"], "columns": lowcov_table["headers"], "style": "Table Style Light 1"})


""" Coverage sheet"""
worksheet_cov.set_column(1, 2, 10)
worksheet_cov.set_column(5, 5, 15)
worksheet_cov.write(0, 0, "Average Coverage per Exon", format_heading)
worksheet_cov.write_row(1, 0, empty_list, format_line)
worksheet_cov.write(2, 0, "Sample: " + str(sample))
worksheet_cov.write(3, 0, "Average coverage of each region in exon-bedfile")

table_area = "A6:G" + str(len(regionscov_table["data"]) + 6)

worksheet_cov.add_table(table_area, {"data": regionscov_table["data"], "columns": regionscov_table["headers"], "style": "Table Style Light 1"})


workbook.close()
