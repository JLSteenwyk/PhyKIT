import pytest

from phykit.errors import PhykitUserError
from phykit.helpers.genomic_coordinates import (
    neighborhoods,
    read_coordinates,
    read_manifest,
    select_rows,
    summarize,
)


def test_manifest_relative_paths_and_duplicate_ids(tmp_path):
    (tmp_path / "a.tre").write_text("(a,b,c,d);")
    manifest = tmp_path / "genes.tsv"
    manifest.write_text("gene_id\ttree\na\ta.tre\n")
    assert read_manifest(manifest) == {"a": str(tmp_path / "a.tre")}
    for body in ("a\ta.tre\na\ta.tre\n", "a\ta.tre\nb\ta.tre\n",
                 "a\tabsent\n", "", "a\n", "\ta.tre\n", "a\ta.tre\textra\n"):
        manifest.write_text("gene_id\ttree\n" + body)
        with pytest.raises(PhykitUserError):
            read_manifest(manifest)


def test_coordinate_merge_and_separate_references(tmp_path):
    bed = tmp_path / "genes.bed"
    bed.write_text("# note\ntrack name=genes\nchr1\t10\t20\ta\nchr1\t12\t30\ta\n")
    rows, diagnostics = read_coordinates([("r1", "bed", bed), ("r2", "bed", bed)])
    assert len(rows) == 2
    assert rows[0]["start"] == 10 and rows[0]["end"] == 30
    assert rows[0]["anchor"] == 19
    assert diagnostics["merged_records"] == {"r1": 1, "r2": 1}
    bed.write_text("chr1\t10\t20\ta\nchr2\t10\t20\ta\n")
    with pytest.raises(PhykitUserError):
        read_coordinates([("r", "bed", bed)])


def test_coordinate_tsv(tmp_path):
    path = tmp_path / "genes.tsv"
    path.write_text("gene_id\tchromosome\tstart\tend\na\tchr1\t0\t1\n")
    rows, _ = read_coordinates([("r", "tsv", path)])
    assert rows[0]["anchor"] == 0


@pytest.mark.parametrize("body", [
    "chr1\t-1\t2\ta\n", "chr1\t2\t2\ta\n", "chr1\tx\t2\ta\n",
    "chr1\t0\t2\n", "chr1\t0\t2\t.\n", "\t0\t2\ta\n", "",
])
def test_invalid_bed(tmp_path, body):
    path = tmp_path / "bad.bed"
    path.write_text(body)
    with pytest.raises(PhykitUserError):
        read_coordinates([("r", "bed", path)])


def row(gene, anchor, category="topology_1", reference="r", chromosome="chr1"):
    return {"gene_id": gene, "anchor": anchor, "start": anchor, "end": anchor + 1,
                "classification": category, "reference": reference, "chromosome": chromosome}


def test_physical_windows_empty_boundary_and_denominators():
    genes = [row("a", 9), row("b", 10, "unresolved"), row("c", 35, "topology_2")]
    windows = neighborhoods(genes, window_bp=10)
    assert [w["total_genes"] for w in windows] == [1, 1, 0, 1]
    assert windows[0]["topology_1_proportion"] == 1
    assert windows[1]["topology_1_proportion"] is None
    assert windows[-1]["end"] == 36
    assert sum(w["total_genes"] for w in windows) == 3
    assert summarize(genes + genes)["total_genes"] == 3


def test_gene_windows_ties_and_references():
    genes = [row("a", 1), row("b", 1), row("c", 3), row("a", 1, reference="s")]
    windows = neighborhoods(genes, window_genes=2)
    assert [w["total_genes"] for w in windows] == [2, 1, 1]
    assert [w["reference"] for w in windows] == ["r", "r", "s"]


def test_interval_and_chromosome_selection():
    genes = [row("a", 9), row("b", 10), row("c", 20), row("d", 11, chromosome="chr2")]
    selected = select_rows(genes, ["chr1"], (10, 20))
    assert [r["gene_id"] for r in selected] == ["b"]
    windows = neighborhoods(selected, window_bp=4, interval=(10, 20))
    assert [(w["start"], w["end"]) for w in windows] == [(10, 14), (14, 18), (18, 20)]
    assert [w["total_genes"] for w in windows] == [1, 0, 0]


@pytest.mark.parametrize("kwargs", [{}, {"window_bp": 0}, {"window_genes": -1},
                                    {"window_bp": 5, "window_genes": 2}])
def test_invalid_windows(kwargs):
    with pytest.raises(PhykitUserError):
        neighborhoods([], **kwargs)


def test_bad_headers_files_and_interval(tmp_path):
    path = tmp_path / "bad"
    for content in ("gene_id\ttree\ttree\n", "wrong\theader\n"):
        path.write_text(content)
        with pytest.raises(PhykitUserError):
            read_manifest(path)
    for fmt in ("bed", "tsv", "gff"):
        with pytest.raises(PhykitUserError):
            read_coordinates([("r", fmt, tmp_path / "missing")])
    with pytest.raises(PhykitUserError):
        select_rows([], interval=(3, 1))
    with pytest.raises(PhykitUserError):
        neighborhoods([row("a", 1000001)], window_bp=1)
