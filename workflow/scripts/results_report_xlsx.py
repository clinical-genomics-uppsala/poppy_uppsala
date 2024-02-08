#!/bin/python3

import xlsxwriter
import gzip
from pysam import VariantFile


# VEP fields in list to get index
def index_vep(variantfile):
    csq_index = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csq_index


# Extract table columns from vcf record
# af, dp, gene, transcipt, coding_name, ensp, consequence, cosmic, clinical, rs, max_pop_af, max_pops, filter_flag, svlen
def extract_vcf_values(record):
    return_dict = {}
    csq = record.info["CSQ"][0].split("|")

    try:
        return_dict["af"] = float(record.info["AF"][0])
    except KeyError:
        pass

    try:
        return_dict["dp"] = int(record.info["DP"])
    except KeyError:
        pass

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        pass
    
    return_dict["artifact_callers"] = record.info["Artifact"]
    return_dict["artifact_median"] = record.info["ArtifactMedian"]
    return_dict["artifact_nr_sd"] = record.info["ArtifactNrSD"]
    return_dict["gene"] = csq[csq_index.index("SYMBOL")]
    return_dict["transcript"] = csq[csq_index.index("HGVSc")].split(":")[0]

    if len(csq[csq_index.index("HGVSc")].split(":")) > 1:
        return_dict["coding_name"] = csq[csq_index.index("HGVSc")].split(":")[1]
    else:
        return_dict["coding_name"] = ""
    return_dict["ensp"] = csq[csq_index.index("HGVSp")]
    return_dict["consequence"] = csq[csq_index.index("Consequence")]

    existing = csq[csq_index.index("Existing_variation")].split("&")
    cosmic_list = [cosmic for cosmic in existing if cosmic.startswith("CO")]
    if len(cosmic_list) == 0:
        return_dict["cosmic"] = ""
    else:
        return_dict["cosmic"] = ", ".join(cosmic_list)

    return_dict["clinical"] = csq[csq_index.index("CLIN_SIG")]

    rs_list = [rs for rs in existing if rs.startswith("rs")]
    if len(rs_list) == 0:
        return_dict["rs"] = ""
    else:
        return_dict["rs"] = ", ".join(rs_list)
    return_dict["max_pop_af"] = csq[csq_index.index("MAX_AF")]
    return_dict["max_pops"] = csq[csq_index.index("MAX_AF_POPS")]
    return_dict["filter_flag"] = ",".join(record.filter.keys())
    # normal freq or, background

    return return_dict


""" SNVs table """
vcf_file = VariantFile(snakemake.input.vcf)
sample = list(vcf_file.header.samples)[0]
sequenceid = snakemake.params.sequenceid
csq_index = index_vep(vcf_file)
snv_table = []
for record in vcf_file.fetch():
    record_values = extract_vcf_values(record)
    if record_values["af"] > 0.01:
        outline = [
            sequenceid,
            sample,
            record_values["gene"],
            record.contig,
            int(record.pos),
            record.ref,
            record.alts[0],
            record_values["af"],
            record_values["dp"],
            record_values["transcript"],
            record_values["coding_name"],
            record_values["ensp"],
            record_values["consequence"],
            record_values["cosmic"],
            record_values["clinical"],
            record_values["rs"],
            record_values["max_pop_af"],
            record_values["max_pops"],
            record_values["filter_flag"],
            record_values["artifact_median"],
            record_values["artifact"],
        ]
        snv_table.append(outline)

""" Pindel table """
pindel_genes = ""

pindel_file = VariantFile(snakemake.input.pindel)
sample = list(pindel_file.header.samples)[0]
sequenceid = snakemake.params.sequenceid
csq_index = index_vep(pindel_file)

pindel_table = []
for record in pindel_file.fetch():
    record_values = extract_vcf_values(record)
    if record_values["af"] > 0.01:
        outline = [
            sequenceid,
            sample,
            record_values["gene"],
            record.contig,
            int(record.pos),
            record.ref,
            record.alts[0],
            record_values["svlen"],
            record_values["af"],
            record_values["dp"],
            record_values["transcript"],
            record_values["coding_name"],
            record_values["ensp"],
            record_values["consequence"],
            record_values["cosmic"],
            record_values["clinical"],
            record_values["rs"],
            record_values["max_pop_af"],
            record_values["max_pops"],
            record_values["filter_flag"],
        ]
        pindel_table.append(outline)


""" .xlsx file with sheets """
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
worksheet_snv = workbook.add_worksheet("SNVs")
worksheet_pindel = workbook.add_worksheet("Pindel")


empty_list = ["", "", "", "", "", ""]
format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_line = workbook.add_format({"top": 1})
format_orange = workbook.add_format({"bg_color": "#ffd280"})
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
worksheet_overview.write_url(8, 0, "internal:'SNVs'!A1", string="SNVs identified")
worksheet_overview.write_url(9, 0, "internal:'Pindel'!A1", string="Pindel results")

## Artefact? Normalpanel?
worksheet_overview.write(23, 0, "Twist Myeloid capture panel TE-91316667_hg38 used")
worksheet_overview.write(
    24,
    0,
    "with the pipeline Poppy ("
    + snakemake.params.poppy_version
    + ") and local uppsala implementations ("
    + snakemake.params.uppsala_version
    + ").",
)
worksheet_overview.write(25, 0, "link to poppy and uppsala repo")

# SNVs sheet
for x in vcf_file.header.records:
    if x.key == "VEP":
        vep_line = x.value

worksheet_snv.set_column(4, 4, 10)
worksheet_snv.set_column(10, 12, 10)
worksheet_snv.write("A1", "Variants found", format_heading)
worksheet_snv.write("A3", "Sample: " + str(sample))
worksheet_snv.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_snv.write("A6", "Databases used: " + vep_line)

worksheet_snv.write("A8", "Filters: ", format_orange)
worksheet_snv.write("B9", "DP_200: Soft filter on depth lower than 200X", format_orange)
worksheet_snv.write("B10", "AD_5: Soft filter variants with few observations (AD lower than 5)", format_orange)
worksheet_snv.write("B11", "AF_0.01: Soft filter variants with low vaf (AF lower than 0.01)", format_orange)
worksheet_snv.write(
    "B12", "PopAF_0.02: : Soft filter germline if >2 % in any population from 1000 genomes, ESP or gnomADe", format_orange
)
worksheet_snv.write(
    "B13",
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic ID on the position",
    format_orange,
)

## Instruction about how to remove filtered variants?
table_area = "A18:U" + str(len(snv_table) + 18)
table_header = [
    {"header": "RunID"},
    {"header": "Sample"},
    {"header": "Gene"},
    {"header": "Chr"},
    {"header": "Pos"},
    {"header": "Ref"},
    {"header": "Alt"},
    {"header": "AF"},
    {"header": "DP"},
    {"header": "Transcript"},
    {"header": "Mutation cds"},
    {"header": "ENSP"},
    {"header": "Consequence"},
    {"header": "COSMIC ids on pos"},
    {"header": "Clinical Significance"},
    {"header": "dbSNP"},
    {"header": "Max Pop AF"},
    {"header": "Max Pop"},
    {"header": "Filter"},
    {"header": "Artifact Medians"},
    {"header": "Artifact calls (Mutect, Vardict, TotNormals)"},
]
worksheet_snv.add_table(table_area, {"data": snv_table, "columns": table_header, "style": "Table Style Light 1"})


worksheet_snv.conditional_format(
    "A17:U" + str(len(snv_table) + 17), {"type": "formula", "criteria": '=LEFT($S16, 4)<>"PASS"', "format": format_orange}
)
# Pindel sheet
worksheet_pindel.set_column(4, 4, 10)
worksheet_pindel.set_column(10, 12, 10)
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
worksheet_pindel.write("B" + str(i + 4), "AF_0.01: Soft filter variants with low vaf (AF lower than 0.01)", format_orange)
worksheet_pindel.write(
    "B" + str(i + 5),
    "PopAF_0.02: : Soft filter germline if >2 % in any population from 1000 genomes, ESP or gnomADe",
    format_orange,
)
worksheet_pindel.write(
    "B" + str(i + 6),
    "Intron: Soft filter intronic variants except if also splice, in cosmic, or in GATA2 or TERT genes or has any Cosmic ID on the position",
    format_orange,
)
i = i + 7

table_area = "A" + str(i) + ":T" + str(len(snv_table) + i)
table_header = [
    {"header": "RunID"},
    {"header": "Sample"},
    {"header": "Gene"},
    {"header": "Chr"},
    {"header": "Pos"},
    {"header": "Ref"},
    {"header": "Alt"},
    {"header": "SV length"},
    {"header": "AF"},
    {"header": "DP"},
    {"header": "Transcript"},
    {"header": "Mutation cds"},
    {"header": "ENSP"},
    {"header": "Consequence"},
    {"header": "COSMIC ids on pos"},
    {"header": "Clinical Significance"},
    {"header": "dbSNP"},
    {"header": "Max Pop AF"},
    {"header": "Max Pop"},
    {"header": "Filter"},
]
worksheet_pindel.add_table(table_area, {"data": pindel_table, "columns": table_header, "style": "Table Style Light 1"})

worksheet_pindel.conditional_format(
    "A17:S" + str(len(pindel_table) + i), {"type": "formula", "criteria": '=LEFT($T16, 4)<>"PASS"', "format": format_orange}
)

# coverage
workbook.close()
