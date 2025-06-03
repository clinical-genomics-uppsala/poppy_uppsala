#!/bin/python3

from results_report_create_tables import *
from datetime import date
from operator import itemgetter
import subprocess
import yaml
import logging

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


def convert_columns_to_letter(nr_columns):
    # Function to convert number of columns to alphabetical coordinates for xlsx-sheets
    if nr_columns < 27:
        letter = chr(nr_columns + 64)
    elif nr_columns < 703:
        i = int((nr_columns - 1) / 26)
        letter = chr(i + 64) + chr(nr_columns - (i * 26) + 64)
    else:
        logging.error(f"Nr columns has to be less than 703, does not support three letter column-index for tables {nr_columns=}")
        sys.exit()
    return letter


vcf = snakemake.input.vcf
pindel = snakemake.input.pindel
sequenceid = snakemake.params.sequenceid
poppy_version = snakemake.params.poppy_version
uppsala_version = snakemake.params.uppsala_version["poppy_uppsala"]["version"]

sample = snakemake.params.sample
sample_type = snakemake.params.sample_type
# Panels dict contain all info for sheets with subsets all SNVs, only need to change here if adding/removing subset sheets
panels = {
    "cll": {"bedfile": snakemake.input.cllbed, "vcf": snakemake.input.cll_vcf},
    "myeloid": {"bedfile": snakemake.input.myeloidbed, "vcf": snakemake.input.myeloid_vcf},
    "hotspot": {"bedfile": snakemake.input.hotspotbed, "vcf": snakemake.input.hotspot_vcf},
}

non_coding_regions = {
    "TERC": ["chr3", 169764300, 169766000, "entire gene + promotor"],
    "GATA2": ["chr3", 128481912, 128483849, "Intron4 (b/w exon4-5) in NM_001145662"],
    "ANKRD26": ["chr10", 27100074, 27100510, "Promotor and exon1"],
    "TP53": ["chr17", 7687366, 7687500, "Exon1 non-coding"],
    "NOTCH1": ["chr9", 136494400, 136496072, "3UTR"],
}
intron_coordinates = {}
for gene in non_coding_regions:
    chrom = non_coding_regions[gene][0]
    if chrom in intron_coordinates:  # If chrom in dict already
        intron_coordinates[chrom].append(non_coding_regions[gene][1:])
    else:
        intron_coordinates[chrom] = [non_coding_regions[gene][1:]]


synonymous_positions = {
    "NM_032638.5(GATA2):c.1416G>A": ["chr3", 128481046, "C", "T"],
    "NM_032638.5(GATA2):c.1023C>T": ["chr3", 128481939, "G", "A"],
    "NM_032638.5(GATA2):c.981G>A": ["chr3", 128483896, "C", "T"],
    "NM_032638.5(GATA2):c.649C>T": ["chr3", 128485949, "G", "A"],
    "NM_032638.5(GATA2):c.351C>G": ["chr3", 128486247, "G", "C"],
    "NM_000546.6(TP53):c.375G>A": ["chr17", 7675994, "C", "T"],
    "NM_000546.6(TP53):c.375G>T": ["chr17", 7675994, "C", "A"],
    "NM_000546.6(TP53):c.375G>C": ["chr17", 7675994, "C", "G"],
    "NM_000546.6(TP53):c.672G>A": ["chr17", 7674859, "C", "T"],
    "NM_000546.6(TP53):c.993G>A": ["chr17", 7673535, "C", "T"],
}


wanted_transcripts = []
with open(snakemake.input.wanted_transcripts) as wanted_file:
    for line in wanted_file:
        wanted_transcripts.append(line.split()[1].split(".")[0])

thresholds = [int(x) for x in snakemake.params.thresholds.split(",")]

for x in VariantFile(vcf).header.records:
    if x.key == "VEP":
        vep_line = x.value

# Add filters from yaml-files
filters_snv = []
filters_pindel = []

for filter_file in [snakemake.params.filter_somatic_hard, snakemake.params.filter_somatic, snakemake.params.filter_pindel]:
    with open(filter_file, "r") as file:
        filters_dict = yaml.safe_load(file)

    for key, items in filters_dict["filters"].items():
        if "hard_filter_somatic" in filter_file:
            filters_snv.append(key + ": " + items["description"])
            filters_pindel.append(key + ": " + items["description"])
        elif "soft_filter_somatic" in filter_file:
            filters_snv.append(items["soft_filter_flag"] + ": " + items["description"])
        elif "soft_filter_pindel" in filter_file:
            filters_pindel.append(items["soft_filter_flag"] + ": " + items["description"])
        else:
            logging.error(f"Unknown filter file, have the filename changed? {filter_file=}")


""" Create data tables to populate excel sheets with """
logging.info(f"Creating data tables to populate excel sheets with")
snv_table = create_snv_table(vcf, sequenceid)
pindel_table = create_pindel_table(pindel, sequenceid)
known_table, known_percent = create_known_variants_table(vcf, pindel, sequenceid)

for panel in panels.keys():
    logging.debug(f"Creating {panel} table")
    panels[panel]["table"] = create_snv_table(panels[panel]["vcf"], sequenceid)
logging.debug(f"{panels=}")


logging.debug(f"Intron and synonymous tables")
intron_table = []
synonymous_table = []
for record in snv_table["data"]:
    if record[4] in intron_coordinates:
        for pair in intron_coordinates[record[4]]:
            if record[5] >= pair[0] and record[5] <= pair[1]:
                intron_table.append(record)

    for position_list in [coordinate for coordinate in synonymous_positions.values() if record[4] == coordinate[0]]:
        if position_list[1] == record[5] and position_list[2] == record[6] and position_list[3] == record[7]:
            synonymous_table.append(record)  # AF filter?


logging.debug(f"List genes in pindel bed")
pindel_genes = []
with open(snakemake.input.pindelbed, "r") as pindel_file:
    for lline in pindel_file:
        line = lline.strip().split("\t")
        pindel_genes.append(line[3])
pindel_genes = list(dict.fromkeys(pindel_genes))

logging.debug(f"Avg cov per exon/coding region from mosdepth {snakemake.input.mosdepth_regions=}")
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

logging.debug(f"Low coverage data {snakemake.input.mosdepth_perbase=}")
lowcov_lines = []
with open(snakemake.input.mosdepth_perbase, "r") as mosdepth_perbase:
    for lline in mosdepth_perbase:
        line = lline.strip().split("\t")
        if int(line[3]) <= int(thresholds[0]) and line[0:4] not in lowcov_lines:
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

for line in lowcov_lines:
    line[3] = int(line[3])
    exons = []
    for bed_line in bed_table:
        # get all exons that cover that low cov line
        if line[0] == bed_line[0] and int(line[1]) >= int(bed_line[1]) and int(line[2]) <= int(bed_line[2]):
            exons.append(bed_line[3])

    # All transcipts in wanted_transcripts file
    if len(exons) > 0:
        wanted_found = []
        for region in exons:
            exon = "_".join(region.split("_")[1:3])
            if exon in wanted_transcripts:
                wanted_found.append(region)
        lowcov_table["data"].append(line + [";".join(wanted_found)] + [";".join(exons)])


logging.debug(f"Overview sheets qc-table values extracted")
coverage = {}
with open(snakemake.input.mosdepth_summary, "r") as summary_file:
    for lline in summary_file:
        line = lline.strip().split("\t")
        if line[0] == "total_region":
            coverage["avg_cov"] = line[3]
        elif line[0] == "chrX_region":
            coverage["chrX_cov"] = line[3]
        elif line[0] == "chrY_region":
            coverage["chrY_cov"] = line[3]

cmd_duplication = "grep -A1 PERCENT_DUPLICATION  " + snakemake.input.picard_dupl + "  | tail -1 | cut -f9"
duplication_rate = (
    float(subprocess.run(cmd_duplication, stdout=subprocess.PIPE, shell="TRUE").stdout.decode("utf-8").strip()) * 100
)

total_breadth = [0, 0, 0]
total_length = 0
with gzip.open(snakemake.input.mosdepth_thresholds, "rt") as threshold_file:
    next(threshold_file)
    for lline in threshold_file:
        line = lline.strip().split("\t")

        length = int(line[2]) - int(line[1])
        total_length += length
        total_breadth[0] += int(line[4])
        total_breadth[1] += int(line[5])
        total_breadth[2] += int(line[6])

thresholds_results = [x / total_length for x in total_breadth]


""" xlsx file with sheets """
logging.info(f"Creating xlsx file {snakemake.output.xlsx}")
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)
worksheet_overview = workbook.add_worksheet("Overview")
if sample.lower() == "hd829":
    worksheet_known = workbook.add_worksheet("Known variants")
else:
    for panel in panels.keys():
        panels[panel]["sheet"] = workbook.add_worksheet(panel.upper())
worksheet_snv = workbook.add_worksheet("SNVs")
worksheet_pindel = workbook.add_worksheet("Pindel")
worksheet_intron = workbook.add_worksheet("Intron")
worksheet_syno = workbook.add_worksheet("Synonymous")
worksheet_lowcov = workbook.add_worksheet("Low Coverage")
worksheet_cov = workbook.add_worksheet("Coverage")
worksheet_qci = workbook.add_worksheet("QCI")


empty_list = ["", "", "", "", "", ""]
format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_line = workbook.add_format({"top": 1})
format_bold = workbook.add_format({"bold": True})
format_orange = workbook.add_format({"bg_color": "#ffd280"})
format_red = workbook.add_format({"font_color": "red"})
format_table_heading = workbook.add_format({"bold": True, "text_wrap": True})

""" Overview """
logging.debug(f"Overview sheet")
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
else:
    for panel in panels.keys():
        worksheet_overview.write_url(i, 0, "internal: '" + panel.upper() + "'!A1", string=panel.upper() + " panel variants")
        i += 1
worksheet_overview.write_url(i, 0, "internal:'SNVs'!A1", string="SNVs identified")
worksheet_overview.write_url(i + 1, 0, "internal:'Pindel'!A1", string="Pindel results")
worksheet_overview.write_url(i + 2, 0, "internal: 'Intron'!A1", string="Intron and non-coding variants")
worksheet_overview.write_url(i + 3, 0, "internal: 'Synonymous'!A1", string="Synonymous variants")
worksheet_overview.write_url(i + 4, 0, "internal:'Low Coverage'!A1", string="Low Coverage regions")
worksheet_overview.write_url(i + 5, 0, "internal:'Coverage'!A1", string="Coverage")
worksheet_overview.write_url(i + 6, 0, "internal: 'QCI'!A1", string="QCI")
i += 8

if sample.lower() == "hd829":
    worksheet_overview.write(i, 0, "Percent of known variants identified:")
    if known_percent < 1:
        worksheet_overview.write(i + 1, 0, str(known_percent * 100) + " %", format_red)
    else:
        worksheet_overview.write(i + 1, 0, str(known_percent * 100) + " %")
    i += 3

worksheet_overview.write_row(
    i,
    0,
    [
        "RunID",
        "DNAnr",
        "Avg. coverage [x]",
        "Duplicationlevel [%]",
        str(thresholds[0]) + "x",
        str(thresholds[1]) + "x",
        str(thresholds[2]) + "x",
    ],
    format_table_heading,
)
worksheet_overview.write_row(i + 1, 0, [sequenceid, sample, coverage["avg_cov"], str(duplication_rate)] + thresholds_results)
i += 3

worksheet_overview.write(i, 0, "Average coverage of regions in 'coding exon' bedfile")
worksheet_overview.write_row(i + 1, 0, ["chrX", coverage["chrX_cov"]])
worksheet_overview.write_row(i + 2, 0, ["chrY", coverage["chrY_cov"]])
i += 4
if sample.lower() != "hd829":
    worksheet_overview.write_url(
        i, 0, "external:" + sample + "_" + sample_type + "_" + sequenceid + "_bamsnap/index.html", string="SNV screenshots"
    )
    i += 2

worksheet_overview.write(i, 0, "Twist Myeloid capture panel TE-96463385_hg38 used")
worksheet_overview.write(
    i + 1,
    0,
    "with the pipeline Poppy (v" + poppy_version + ") and local uppsala implementations (v" + uppsala_version + ").",
)
worksheet_overview.write_url(i + 3, 0, "https://gms-poppy.readthedocs.io/en/latest/", string="Poppy documentation")
worksheet_overview.write_url(
    i + 4, 0, "https://github.com/clinical-genomics-uppsala/poppy_uppsala/tree/develop", string="Uppsala configurations"
)
worksheet_overview.write(i + 6, 0, "Specific program versions can be found in MultiQC report")
i += 8

worksheet_overview.write(i, 0, "Full design bedfile: " + snakemake.input.bedfile)
worksheet_overview.write(i + 1, 0, "Coding exons bedfile: " + snakemake.input.exonbed)
worksheet_overview.write(i + 2, 0, "Artifact panel used: " + snakemake.params.artifact)
worksheet_overview.write(i + 3, 0, "Background panel used for snvs: " + snakemake.params.background)
worksheet_overview.write(i + 4, 0, "Pindel bedfile used: " + snakemake.input.pindelbed)
worksheet_overview.write(i + 5, 0, "Pindel artifact panel used: " + snakemake.params.artifact_pindel)
i += 6
if sample.lower() != "hd829":
    for panel in panels.keys():
        worksheet_overview.write(i, 0, panel.upper() + "bedfile used: " + panels[panel]["bedfile"])
        i += 1


""" Known sheet """
logging.debug(f"Known variants sheet if hd829, else panel sheets {sample.lower()=}")
if sample.lower() == "hd829":
    worksheet_known.set_column(3, 4, 10)
    worksheet_known.set_column(6, 6, 10)
    worksheet_known.write("A1", "Variants known for HD829", format_heading)
    worksheet_known.write("A3", "Sample: " + str(sample))
    i = 5
    column_end = ":" + convert_columns_to_letter(len(known_table["headers"]))
    table_area = "A" + str(i) + column_end + str(len(known_table["data"]) + i)
    worksheet_known.add_table(
        table_area, {"data": known_table["data"], "columns": known_table["headers"], "style": "Table Style Light 1"}
    )
else:
    """CLL and myeloid tables (or other panels)"""
    logging.debug(f"panel sheets {panels.keys()=}")
    for panel in panels.keys():
        logging.debug(f"Create {panel} sheet. {panels[panel]=}")
        worksheet_panel = panels[panel]["sheet"]
        worksheet_panel.set_column(2, 2, 10)
        worksheet_panel.set_column(5, 5, 10)
        worksheet_panel.set_column(11, 13, 10)
        worksheet_panel.write("A1", "Variants found in " + panel.upper() + " gene panel", format_heading)
        worksheet_panel.write("A3", "Sample: " + str(sample))
        worksheet_panel.write("A4", "Reference used: " + str(snakemake.params.ref))
        worksheet_panel.write("A5", "Bedfile used: " + str(panels[panel]["bedfile"]))
        worksheet_panel.write("A7", "Databases used: " + vep_line)

        worksheet_panel.write("A9", "Filters: ", format_orange)
        for i, filter_txt in enumerate(filters_snv):
            i += 10
            worksheet_panel.write("B" + str(i), filter_txt, format_orange)

        i += 2
        worksheet_panel.write_rich_string(
            "A" + str(i),
            "Only variants with ",
            format_bold,
            "> 2 % AF",
            " and filter-flag ",
            format_bold,
            "PASS",
            " shown by default.",
        )
        worksheet_panel.write(
            "A" + str(i + 1),
            "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
            + "You can then use the drop-downs in the header row to filter to your liking.",
        )

        worksheet_panel.write_url(
            "A" + str(i + 3),
            "external:" + sample + "_" + sample_type + "_" + sequenceid + "_bamsnap/index.html",
            string="SNV screenshots",
        )
        worksheet_panel.write("A" + str(i + 4), "Only variants with AF >= 5% and PASS have automated screenshots.")

        i += 6
        column_end = ":" + convert_columns_to_letter(len(panels[panel]["table"]["headers"]))
        if len(panels[panel]["table"]["data"]) > 0:
            table_area = "A" + str(i) + column_end + str(len(panels[panel]["table"]["data"]) + i)
            table_area_data = "A" + str(i + 1) + column_end + str(len(panels[panel]["table"]["data"]) + i)
        else:
            table_area = "A" + str(i) + column_end + str(i + 1)
            table_area_data = "A" + str(i + 1) + column_end + str(i + 1)

        worksheet_panel.add_table(table_area, {"columns": panels[panel]["table"]["headers"], "style": "Table Style Light 1"})
        cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
        worksheet_panel.conditional_format(
            table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange}
        )

        worksheet_panel.autofilter(table_area)
        worksheet_panel.filter_column("A", "Filter != PASS")
        worksheet_panel.filter_column("I", "AF >= 0.02")
        for row_data in panels[panel]["table"]["data"]:
            if row_data[0] == "PASS" and float(row_data[8]) >= 0.02:
                pass
            else:
                worksheet_panel.set_row(i, options={"hidden": True})
            worksheet_panel.write_row(i, 0, row_data)
            i += 1

""" SNVs sheet """
logging.debug(f"SNVs sheet")
worksheet_snv.set_column(2, 2, 10)
worksheet_snv.set_column(5, 5, 10)
worksheet_snv.set_column(11, 13, 10)
worksheet_snv.write("A1", "Variants found", format_heading)
worksheet_snv.write("A3", "Sample: " + str(sample))
worksheet_snv.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_snv.write("A6", "Databases used: " + vep_line)

worksheet_snv.write("A8", "Filters: ", format_orange)
for i, filter_txt in enumerate(filters_snv):
    i += 9
    worksheet_snv.write("B" + str(i), filter_txt, format_orange)

i += 2
worksheet_snv.write_rich_string(
    "A" + str(i), "Only variants with ", format_bold, "> 2 % AF", " and filter-flag ", format_bold, "PASS", " shown by default."
)
worksheet_snv.write(
    "A" + str(i + 1),
    "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
    + "You can then use the drop-downs in the header row to filter to your liking.",
)

worksheet_snv.write_url(
    "A" + str(i + 3),
    "external:" + sample + "_" + sample_type + "_" + sequenceid + "_bamsnap/index.html",
    string="SNV screenshots",
)
worksheet_snv.write("A" + str(i + 4), "Only variants with AF >= 5% and PASS have automated screenshots.")

i += 6
column_end = ":" + convert_columns_to_letter(len(snv_table["headers"]))
if len(snv_table["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(snv_table["data"]) + i)
    table_area_data = "A" + str(i + 1) + column_end + str(len(snv_table["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)
    table_area_data = "A" + str(i + 1) + column_end + str(i + 1)

worksheet_snv.add_table(table_area, {"columns": snv_table["headers"], "style": "Table Style Light 1"})
cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
worksheet_snv.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})


worksheet_snv.autofilter(table_area)
worksheet_snv.filter_column("A", "Filter != PASS")
worksheet_snv.filter_column("I", "AF >= 0.02")
for row_data in snv_table["data"]:
    if row_data[0] == "PASS" and float(row_data[8]) >= 0.02:
        pass
    else:
        worksheet_snv.set_row(i, options={"hidden": True})
    worksheet_snv.write_row(i, 0, row_data)
    i += 1


""" Pindel sheet """
logging.debug(f"Pindel sheet")
worksheet_pindel.set_column(2, 2, 10)
worksheet_pindel.set_column(5, 5, 10)
worksheet_pindel.set_column(11, 13, 10)
worksheet_pindel.write("A1", "Variants found", format_heading)
worksheet_pindel.write("A3", "Sample: " + str(sample))
worksheet_pindel.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_pindel.write("A6", "To limit runtime pindel were used with a specific designfile: " + snakemake.input.pindelbed)
worksheet_pindel.write("A7", "Which includes the following genes: ")
i = 8
for gene in pindel_genes:
    worksheet_pindel.write("C" + str(i), gene)
    i += 1

worksheet_pindel.write("A" + str(i + 1), "Filters: ", format_orange)
for j, filter_txt in enumerate(filters_pindel):
    j += i + 1
    worksheet_pindel.write("B" + str(j), filter_txt, format_orange)
i += 2 + len(filters_pindel)

worksheet_pindel.write_rich_string("A" + str(i), "Only variants with filter-flag ", format_bold, "PASS", " shown by default.")
worksheet_pindel.write(
    "A" + str(i + 1),
    "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
    + "You can then use the drop-downs in the header row to filter to your liking.",
)
i += 3


column_end = ":" + convert_columns_to_letter(len(pindel_table["headers"]))
if len(pindel_table["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(pindel_table["data"]) + i)
    table_area_data = "A" + str(i + 1) + column_end + str(len(pindel_table["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)
    table_area_data = "A" + str(i + 1) + column_end + str(i + 1)

worksheet_pindel.add_table(table_area, {"columns": pindel_table["headers"], "style": "Table Style Light 1"})
cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
worksheet_pindel.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

worksheet_pindel.autofilter(table_area)
worksheet_pindel.filter_column("A", "Filter != PASS")
for row_data in pindel_table["data"]:
    if row_data[0] == "PASS":
        pass
    else:
        worksheet_pindel.set_row(i, options={"hidden": True})
    worksheet_pindel.write_row(i, 0, row_data)
    i += 1


""" Intron variants """
logging.debug(f"Intron sheet")
worksheet_intron.set_column(2, 2, 10)
worksheet_intron.set_column(3, 5, 10)
worksheet_intron.set_column(11, 13, 10)

worksheet_intron.write("A1", "Intron and non-coding variants in selected regions", format_heading)
worksheet_intron.write("A3", "Sample: " + str(sample))
worksheet_intron.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_intron.write("A6", "Intron variants found in the following regions only: ")
i = 6
for gene in non_coding_regions:
    worksheet_intron.write_row(i, 1, [gene] + non_coding_regions[gene])
    i += 1

i += 3
column_end = ":" + convert_columns_to_letter(len(snv_table["headers"]))
table_area = "A" + str(i) + column_end + str(len(intron_table) + i)
worksheet_intron.add_table(table_area, {"data": intron_table, "columns": snv_table["headers"], "style": "Table Style Light 1"})


""" Synonymous variants """
logging.debug(f"Synonymous sheet")
worksheet_syno.write("A1", "Synonymous variants at selected positions", format_heading)
worksheet_syno.write("A3", "Sample: " + str(sample))
worksheet_syno.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_syno.write("A6", "Synonymous variants found in the following positions and matching ALT only: ")
i = 6
for c_name, values in synonymous_positions.items():
    worksheet_syno.write_row(i, 1, [c_name] + values)
    i += 1

i += 3
if synonymous_table == 0:
    synonymous_table.append([""] * 24)
column_end = ":" + convert_columns_to_letter(len(snv_table["headers"]))
table_area = "A" + str(i) + column_end + str(len(synonymous_table) + i)
worksheet_syno.add_table(table_area, {"data": synonymous_table, "columns": snv_table["headers"], "style": "Table Style Light 1"})


""" Low coverage """
logging.debug(f"Low coverage sheet")
worksheet_lowcov.set_column(1, 2, 10)
worksheet_lowcov.set_column(4, 5, 25)

worksheet_lowcov.write(0, 0, "Mosdepth low coverage analysis", format_heading)
worksheet_lowcov.write_row(1, 0, empty_list, format_line)
worksheet_lowcov.write(2, 0, "Sample: " + str(sample))
worksheet_lowcov.write(3, 0, "Gene regions with coverage lower than " + str(thresholds[0]) + "x.")

column_end = ":" + convert_columns_to_letter(len(lowcov_table["headers"]))
table_area = "A6" + column_end + str(len(lowcov_table["data"]) + 6)
worksheet_lowcov.add_table(
    table_area, {"data": lowcov_table["data"], "columns": lowcov_table["headers"], "style": "Table Style Light 1"}
)


""" Coverage sheet"""
logging.debug(f"Coverage sheet")
worksheet_cov.set_column(1, 2, 10)
worksheet_cov.set_column(5, 5, 15)
worksheet_cov.write(0, 0, "Average Coverage per Exon", format_heading)
worksheet_cov.write_row(1, 0, empty_list, format_line)
worksheet_cov.write(2, 0, "Sample: " + str(sample))
worksheet_cov.write(3, 0, "Average coverage of each region in exon-bedfile")

column_end = ":" + convert_columns_to_letter(len(regionscov_table["headers"]))
table_area = "A6" + column_end + str(len(regionscov_table["data"]) + 6)

worksheet_cov.add_table(
    table_area, {"data": regionscov_table["data"], "columns": regionscov_table["headers"], "style": "Table Style Light 1"}
)


""" QCI sheet """
logging.debug(f"QCI sheet")
qci_table_header = [
    "DNA nr",
    "Chromosome",
    "Position",
    "Gene Region",
    "Gene Symbol",
    "Transcript ID",
    "Transcript Variant",
    "Protein Variant",
    "Variant Findings",
    "Sample Genotype Quality",
    "Read Depth",
    "Allele Fraction",
    "Translation Impact",
    "dbSNP ID",
    "1000 Genomes Frequency",
    "ExAC Frequency",
    "HGMD",
    "COSMIC ID",
    "Artefacts_without_ASXL1",
    "ASXL1_variant_filter",
]
worksheet_qci.set_column("C:C", 10)
worksheet_qci.write("A1", "Results from QCI", format_heading)
worksheet_qci.write_row("A2", empty_list, format_line)

worksheet_qci.write("A5", "Analysen utfördes i enlighet med dokumentationen.")
worksheet_qci.write("A6", "Eventuella avvikelser: ")

worksheet_qci.write_row(9, 0, qci_table_header, format_table_heading)

workbook.close()
logging.info(f"All done!")
