"""
Unit tests for quartet_pie (quartet concordance pie chart visualization).

Tests initialization, process_args, plot creation, and JSON output.
Ground truth gCF/gDF values validated against manual bipartition
matching computation on sample data.
"""
import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.quartet_pie import QuartetPie


@pytest.fixture
def args(tmp_path):
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        gene_trees="tests/sample_files/gene_trees_simple.nwk",
        output=str(tmp_path / "qpie.png"),
        annotate=False,
    )


class TestQuartetPieInit:
    def test_init_sets_fields(self, args):
        qp = QuartetPie(args)
        assert qp.tree_file_path == args.tree
        assert qp.gene_trees_path == args.gene_trees
        assert qp.output_path == args.output
        assert qp.annotate is False
        assert qp.json_output is False

    def test_process_args_no_gene_trees(self, tmp_path):
        args = Namespace(
            tree="t.tre", output=str(tmp_path / "out.png"),
        )
        qp = QuartetPie(args)
        assert qp.gene_trees_path is None


class TestQuartetPiePlot:
    def test_creates_png(self, args):
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_creates_pdf(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            output=str(tmp_path / "qpie.pdf"),
            annotate=False,
        )
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()

    def test_annotate_flag(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            output=str(tmp_path / "qpie_annot.png"),
            annotate=True,
        )
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()


class TestQuartetPieJson:
    def test_json_output(self, mocker, args):
        args.json = True
        qp = QuartetPie(args)
        mocked_json = mocker.patch("phykit.services.tree.quartet_pie.print_json")
        qp.run()

        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_gene_trees"] == 10
        assert payload["input_mode"] == "native"
        assert len(payload["nodes"]) == 5

        # Verify ground truth for bear+raccoon node
        br_node = next(
            n for n in payload["nodes"]
            if set(n["node_tips"]) == {"bear", "raccoon"}
        )
        assert br_node["gCF"] == 1.0
        assert br_node["concordant_count"] == 7


class TestQuartetPieValidation:
    def test_too_few_tips_exits(self, tmp_path):
        from Bio import Phylo
        from io import StringIO

        # Write a 3-tip tree
        tree_path = tmp_path / "small.tre"
        tree_path.write_text("(A:1,B:1,C:1);")

        args = Namespace(
            tree=str(tree_path),
            gene_trees=None,
            output=str(tmp_path / "out.png"),
            annotate=False,
        )
        qp = QuartetPie(args)
        with pytest.raises(SystemExit):
            qp.run()
