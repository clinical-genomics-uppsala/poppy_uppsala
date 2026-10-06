__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2026, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

from pysam import VariantFile, VariantHeader

from cnv_large_calls_table import get_cytoband, get_large_cnvs, read_cytobands

CYTOBANDS = {
    "chr1": [(0, 2300000, "p36.33"), (2300000, 5300000, "p36.32"), (5300000, 7100000, "p36.31")],
}


def test_get_cytoband():
    assert get_cytoband(CYTOBANDS, "chr1", 100, 2000) == "1p36.33"
    assert get_cytoband(CYTOBANDS, "chr1", 2300001, 2400000) == "1p36.32"
    assert get_cytoband(CYTOBANDS, "chr1", 2300000, 5300001) == "1p36.33-p36.31"
    assert get_cytoband(CYTOBANDS, "chr2", 100, 2000) == ""


def test_read_cytobands(tmp_path):
    cytoband_file = tmp_path / "cytoBand.txt"
    cytoband_file.write_text("chr1\t2300000\t5300000\tp36.32\tgpos25\nchr1\t0\t2300000\tp36.33\tgneg\n")
    assert read_cytobands(cytoband_file) == {"chr1": [(0, 2300000, "p36.33"), (2300000, 5300000, "p36.32")]}


def test_get_large_cnvs(tmp_path):
    header = VariantHeader()
    header.contigs.add("chr1", length=7100000)
    header.add_meta("INFO", items=[("ID", "END"), ("Number", 1), ("Type", "Integer"), ("Description", "End")])
    header.add_meta("INFO", items=[("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "Length")])
    header.add_meta("INFO", items=[("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "Type")])
    header.add_meta("INFO", items=[("ID", "CALLER"), ("Number", 1), ("Type", "String"), ("Description", "Caller")])
    header.add_meta("INFO", items=[("ID", "CORR_CN"), ("Number", 1), ("Type", "Float"), ("Description", "CN")])
    header.add_meta("INFO", items=[("ID", "Normal_AF"), ("Number", 1), ("Type", "Float"), ("Description", "AF")])
    vcf_file = tmp_path / "cnv.vcf"
    calls = [
        # (start, end, svtype, normal_af)
        (1000, 2500000, "DEL", None),  # kept, no Normal_AF annotation
        (3000000, 3050000, "DUP", 0.0),  # too short
        (3000000, 5000000, "COPY_NORMAL", 0.0),  # copy neutral
        (5400000, 7000000, "DUP", 0.5),  # artifact
        (5400000, 7000000, "DUP", 0.1),  # kept
    ]
    with VariantFile(str(vcf_file), "w", header=header) as vcf:
        for start, end, svtype, normal_af in calls:
            info = {"SVTYPE": svtype, "SVLEN": end - start, "CALLER": "cnvkit", "CORR_CN": 1.234}
            if normal_af is not None:
                info["Normal_AF"] = normal_af
            vcf.write(vcf.new_record(contig="chr1", start=start - 1, stop=end, alleles=("N", f"<{svtype}>"), info=info))

    rows = get_large_cnvs(str(vcf_file), CYTOBANDS, min_length=100000, max_normal_af=0.15)

    assert [(row["Start"], row["Type"]) for row in rows] == [(1000, "DEL"), (5400000, "DUP")]
    assert rows[0] == {
        "Chromosome": "chr1",
        "Start": 1000,
        "End": 2500000,
        "Cytoband": "1p36.33-p36.32",
        "Type": "DEL",
        "Copy number": 1.23,
        "Length": 2499000,
        "Caller": "cnvkit",
    }
