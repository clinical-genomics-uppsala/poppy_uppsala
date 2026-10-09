__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2026, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

from pysam import VariantFile, VariantHeader

from cnv_large_calls_table import chrom_sort_key, get_cytoband, get_large_cnvs, read_cytobands

CYTOBANDS = {
    "chr1": [(0, 2300000, "p36.33"), (2300000, 5300000, "p36.32"), (5300000, 7100000, "p36.31")],
}


def test_get_cytoband():
    assert get_cytoband(CYTOBANDS, "chr1", 100, 2000) == "1p36.33"
    assert get_cytoband(CYTOBANDS, "chr1", 2300001, 2400000) == "1p36.32"
    assert get_cytoband(CYTOBANDS, "chr1", 2300000, 5300001) == "1p36.33-p36.31"
    assert get_cytoband(CYTOBANDS, "chr2", 100, 2000) == ""


def test_chrom_sort_key():
    chroms = ["chrX", "chr10", "chr2", "chrY", "chr1", "chr22", "chr11"]
    assert sorted(chroms, key=chrom_sort_key) == ["chr1", "chr2", "chr10", "chr11", "chr22", "chrX", "chrY"]


def test_read_cytobands(tmp_path):
    cytoband_file = tmp_path / "cytoBand.txt"
    cytoband_file.write_text("chr1\t2300000\t5300000\tp36.32\tgpos25\nchr1\t0\t2300000\tp36.33\tgneg\n")
    assert read_cytobands(cytoband_file) == {"chr1": [(0, 2300000, "p36.33"), (2300000, 5300000, "p36.32")]}


def write_cnv_vcf(vcf_file, calls):
    """Write calls given as (chrom, start, end, svtype, normal_af, baf) tuples to a CNV VCF."""
    header = VariantHeader()
    for contig in ["chr1", "chr10", "chr2", "chrX"]:
        header.contigs.add(contig, length=7100000)
    header.add_meta("INFO", items=[("ID", "END"), ("Number", 1), ("Type", "Integer"), ("Description", "End")])
    header.add_meta("INFO", items=[("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "Length")])
    header.add_meta("INFO", items=[("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "Type")])
    header.add_meta("INFO", items=[("ID", "CALLER"), ("Number", 1), ("Type", "String"), ("Description", "Caller")])
    header.add_meta("INFO", items=[("ID", "CORR_CN"), ("Number", 1), ("Type", "Float"), ("Description", "CN")])
    header.add_meta("INFO", items=[("ID", "Normal_AF"), ("Number", 1), ("Type", "Float"), ("Description", "AF")])
    header.add_meta("INFO", items=[("ID", "BAF"), ("Number", 1), ("Type", "Float"), ("Description", "BAF")])
    with VariantFile(str(vcf_file), "w", header=header) as vcf:
        for chrom, start, end, svtype, normal_af, baf in calls:
            info = {"SVTYPE": svtype, "SVLEN": end - start, "CALLER": "cnvkit", "CORR_CN": 1.234}
            if normal_af is not None:
                info["Normal_AF"] = normal_af
            if baf is not None:
                info["BAF"] = baf
            vcf.write(vcf.new_record(contig=chrom, start=start - 1, stop=end, alleles=("N", f"<{svtype}>"), info=info))


def test_get_large_cnvs(tmp_path):
    vcf_file = tmp_path / "cnv.vcf"
    write_cnv_vcf(vcf_file, [
        # (chrom, start, end, svtype, normal_af, baf), in lexicographic chromosome order
        ("chr1", 1000, 2500000, "DEL", None, 0.123),  # kept, no Normal_AF annotation
        ("chr1", 3000000, 3050000, "DUP", 0.0, None),  # too short
        ("chr1", 3000000, 5000000, "COPY_NORMAL", 0.0, None),  # copy neutral
        ("chr1", 5400000, 7000000, "DUP", 0.5, None),  # artifact
        ("chr1", 5400000, 7000000, "DUP", 0.1, None),  # kept
        ("chr10", 1000, 2500000, "DEL", 0.0, None),  # kept
        ("chr2", 1000, 2500000, "DEL", 0.0, None),  # kept
        ("chrX", 1000, 2500000, "DUP", 0.0, None),  # kept
    ])

    rows = get_large_cnvs(str(vcf_file), CYTOBANDS, min_length=100000, max_normal_af=0.15)

    assert [(row["Chromosome"], row["Start"]) for row in rows] == [
        ("chr1", 1000),
        ("chr1", 5400000),
        ("chr2", 1000),
        ("chr10", 1000),
        ("chrX", 1000),
    ]
    assert rows[0] == {
        "Chromosome": "chr1",
        "Start": 1000,
        "End": 2500000,
        "Cytoband": "1p36.33-p36.32",
        "Type": "DEL",
        "Copy number": 1.23,
        "BAF": 0.12,
        "Length": 2499000,
        "Caller": "cnvkit",
    }
    assert rows[1]["BAF"] == ""


def test_get_large_cnvs_excludes_small_aberrations(tmp_path):
    vcf_file = tmp_path / "cnv.vcf"
    write_cnv_vcf(vcf_file, [
        # (chrom, start, end, svtype, normal_af, baf), SVLEN = end - start
        ("chr1", 1000, 2000, "DEL", 0.0, None),  # 1 kb, well below min_length
        ("chr1", 3000000, 3099999, "DUP", 0.0, None),  # min_length - 1
        ("chr1", 4000000, 4100000, "DEL", 0.0, None),  # exactly min_length, included
        ("chr1", 5000000, 5100001, "DUP", 0.0, None),  # min_length + 1, kept
    ])

    rows = get_large_cnvs(str(vcf_file), CYTOBANDS, min_length=100000, max_normal_af=0.15)

    assert [(row["Start"], row["Length"]) for row in rows] == [(4000000, 100000),(5000000, 100001)]
