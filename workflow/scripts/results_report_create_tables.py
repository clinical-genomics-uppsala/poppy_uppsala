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


# Extract table columns from vcf records
def extract_vcf_values(record, csq_index):
    return_dict = {}
    csq = record.info["CSQ"][0].split("|")

    try:
        return_dict["af"] = float(record.info["AF"][0])
    except KeyError:
        return_dict["af"] = int(record.samples[0].get("AD")[1]) / sum(record.samples[0].get("AD"))

    try:
        return_dict["dp"] = int(record.info["DP"])
    except KeyError:
        return_dict["dp"] = sum(record.samples[0].get("AD"))

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        pass

    try:
        return_dict["artifact_callers"] = ";".join(record.info["Artifact"])
        return_dict["artifact_median"] = ";".join([str(round(float(x), 3)) for x in record.info["ArtifactMedian"]])
        return_dict["artifact_nr_sd"] = ";".join(record.info["ArtifactNrSD"])
    except KeyError:
        pass

    try:
        return_dict["callers"] = ";".join(record.info["CALLERS"])
    except KeyError:
        pass

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
    if return_dict["af"] >= 0.05:
        return_dict["igv"] = "pathpath"
    else:
        return_dict["igv"] = ""
    # normal freq or, background

    return return_dict


def create_snv_table(vcf_input, sequenceid):
    vcf_file = VariantFile(vcf_input)
    sample = list(vcf_file.header.samples)[0]
    csq_index = index_vep(vcf_file)

    snv_table = {"data": [], "headers": []}
    snv_table["headers"] = [
        {"header": "Filter"},
        {"header": "RunID"},
        {"header": "DNAnr"},
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
        {"header": "Artifact Medians (Mutect; Vardict)"},
        {"header": "Artifact calls (Mutect; Vardict; TotNormals)"},
        {"header": "Callers"},
    ]
    for record in vcf_file.fetch():
        record_values = extract_vcf_values(record, csq_index)
        if record_values["af"] > 0.01:
            outline = [
                record_values["filter_flag"],
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
                record_values["artifact_median"],
                record_values["artifact_callers"],
                record_values["callers"],
            ]
            snv_table["data"].append(outline)
    return snv_table


def create_pindel_table(vcf_input, sequenceid):
    pindel_file = VariantFile(vcf_input)
    sample = list(pindel_file.header.samples)[0]
    csq_index = index_vep(pindel_file)

    pindel_table = {"data": [], "headers": []}
    pindel_table["headers"] = [
        {"header": "Filter"},
        {"header": "RunID"},
        {"header": "DNAnr"},
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
    ]
    for record in pindel_file.fetch():
        record_values = extract_vcf_values(record, csq_index)
        if record_values["af"] > 0.01:
            outline = [
                record_values["filter_flag"],
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
            ]
            pindel_table["data"].append(outline)
    return pindel_table


def create_known_variants_table(vcf_input, pindel_input, sequenceid):
    known_number = 0
    vcf_file = VariantFile(vcf_input)
    pindel_file = VariantFile(pindel_input)
    sample = list(vcf_file.header.samples)[0]
    known = [
        ["ABL1", "T315I", "c.1001C>T", "chr9", 130872896, "C", "T", "COSM12560", "COSV59323790", 0.0523],
        ["ASXL1", "G646Wfs*12", "c.1934dup", "chr20", 32434638, "A", "AG", "COSM34210", "COSV60102155", 0.3966],
        ["ASXL1", "W796C", "c.2388G>T", "chr20", 32435100, "G", "T", "COSM1681610", "COSV60115219", 0.0536],
        ["BCOR", "Q1208Tfs*8", "c.3621dup", "chrX", 40063833, "G", "GT", "COSM1683571", "COSV60698892", 0.6844],
        ["CBL", "S403F", "c.1208C>T", "chr11", 119278278, "C", "T", "COSM1676499", "COSV50649507", 0.0562],
        ["DNMT3A", "R882C", "c.2644C>T", "chr2", 25234374, "G", "A", "COSM53042", "COSV53036332", 0.0626],
        ["EZH2", "R418Q", "c.1253G>A", "chr7", 148817379, "C", "T", "COSM3259655", "COSV57458095", 0.0476],
        ["FLT3", "D835Y", "c.2503G>T", "chr13", 28018505, "C", "A", "COSM783", "COSV54042116", 0.0471],
        ["FLT3", "ITD300", "N/A", "chr13", 28033909, "N/A", "N/A", "N/A", "N/A", 0.05],
        ["GATA1", "Q119*", "c.355C>T", "chrX", 48791978, "C", "T", "N/A", "N/A", 0.0805],
        ["GATA2", "G200Vfs*18", "c.599del", "chr3", 128485998, "AC", "A", "COSM1418772", "COSV62003348", 0.3027],
        ["IDH1", "R132C", "c.394C>T", "chr2", 208248389, "G", "A", "COSM28747", "COSV61615256", 0.0442],
        ["IDH2", "R172K", "c.515G>A", "chr15", 90088606, "C", "T", "COSM33733", "COSV57468734", 0.0262],
        ["JAK2", "F537-K539>L", "c.1611_1616del", "chr9", 5070021, "TTCACAA", "T", "COSM24437", "COSV67579858", 0.0401],
        ["JAK2", "V617F", "c.1849G>T", "chr9", 5073770, "G", "T", "COSM12600", "COSV67569051", 0.0419],
        ["KRAS", "G13D", "c.38G>A", "chr12", 25245347, "C", "T", "COSM532", "COSV55497388", 0.3604],
        ["NPM1", "W288Cfs*12", "c.860_863dup", "chr5", 171410539, "C", "CTCTG", "COSM17559", "COSV51542664", 0.0389],
        ["NRAS", "Q61L", "c.182A>T", "chr1", 114713908, "T", "A", "COSM583", "COSV54736624", 0.094],
        ["RUNX1", "M267I", "c.801G>A", "chr21", 34834414, "C", "T", "COSM1681955", "COSV55866591", 0.299],
        ["SF3B1", "G740E", "c.2219G>A", "chr2", 197401989, "C", "T", "COSM133120", "COSV59205460", 0.0314],
        ["TET2", "R1261H", "c.3782G>A", "chr4", 105243757, "G", "A", "COSM211643", "COSV54396706", 0.0494],
        ["TP53", "S241F", "c.722C>T", "chr17", 7674241, "G", "A", "COSM10812", "COSV52661688", 0.0582],
    ]
    known_table = {"data": [], "headers": []}
    known_table["headers"] = [
        {"header": "RunID"},
        {"header": "DNAnr"},
        {"header": "Gene"},
        {"header": "Variant AA"},
        {"header": "CDS mutation"},
        {"header": "Chr"},
        {"header": "Pos"},
        {"header": "Ref"},
        {"header": "Alt"},
        {"header": "Legacy ID"},
        {"header": "Genomic mutation ID"},
        {"header": "Expected AF"},
        {"header": "AF"},
        {"header": "DP"},
        {"header": "SV Length"},
    ]
    for known_line in known:
        if known_line[1] == "ITD300":
            for record in pindel_file.fetch(known_line[3], known_line[4] - 1, known_line[4]):
                af = int(record.samples[0].get("AD")[1]) / sum(record.samples[0].get("AD"))
                dp = sum(record.samples[0].get("AD"))
                outline = [sequenceid, sample] + known_line + [float(af), int(dp), record.info["SVLEN"]]
                known_table["data"].append(outline)
                known_number += 1
        else:
            for record in vcf_file.fetch(known_line[3], known_line[4] - 1, known_line[4]):
                if record.alts[0] == known_line[6] and record.ref == known_line[5]:
                    outline = [sequenceid, sample] + known_line + [float(record.info["AF"][0]), int(record.info["DP"])]
                    known_table["data"].append(outline)
                    known_number += 1

    return known_table, known_number / len(known)
