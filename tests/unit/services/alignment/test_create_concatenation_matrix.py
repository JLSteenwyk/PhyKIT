from argparse import Namespace
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from phykit.services.alignment.create_concatenation_matrix import CreateConcatenationMatrix
import phykit.services.alignment.create_concatenation_matrix as ccm_module


def _write_fasta(path: Path, records):
    lines = []
    for taxon, seq in records:
        lines.append(f">{taxon}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def args():
    return Namespace(alignment_list="/some/path/to/file", prefix="some_prefix")


class TestCreateConcatenationMatrix:
    def test_init_sets_alignment_list_path(self, args):
        ccm = CreateConcatenationMatrix(args)
        assert ccm.alignment_list_path == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
        assert ccm.json_output is False
        assert ccm.plot_occupancy is False
        assert ccm.plot_output is None
        assert ccm.threshold == 0
        assert ccm.plot_config is not None

    def test_process_args_reads_optional_flags(self):
        parsed = CreateConcatenationMatrix(
            Namespace(
                alignment_list="x.list",
                prefix="x",
                json=True,
                plot_occupancy=True,
                plot_output="occ.png",
            )
        ).process_args(
            Namespace(
                alignment_list="x.list",
                prefix="x",
                json=True,
                plot_occupancy=True,
                plot_output="occ.png",
            )
        )
        assert parsed["json_output"] is True
        assert parsed["plot_occupancy"] is True
        assert parsed["plot_output"] == "occ.png"

    def test_read_alignment_paths_missing_file_exits(self, args):
        ccm = CreateConcatenationMatrix(args)
        with pytest.raises(SystemExit) as exc:
            ccm.read_alignment_paths("/definitely/not/found.txt")
        assert exc.value.code == 2

    def test_get_taxa_from_alignment(self, tmp_path, args):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("B", "A-GT")])
        ccm = CreateConcatenationMatrix(args)
        taxa = ccm._get_taxa_from_alignment(str(fasta))
        assert taxa == {"A", "B"}

    def test_get_taxa_names_sequential(self, tmp_path, args):
        fastas = []
        for idx, entries in enumerate(
            [
                [("A", "AAAA"), ("B", "CCCC")],
                [("A", "GGGG"), ("C", "TTTT")],
            ]
        ):
            fasta = tmp_path / f"gene_{idx}.fa"
            _write_fasta(fasta, entries)
            fastas.append(str(fasta))

        ccm = CreateConcatenationMatrix(args)
        taxa = ccm.get_taxa_names(fastas)
        assert taxa == ["A", "B", "C"]

    def test_get_taxa_names_parallel_fallback(self, monkeypatch, args):
        class FailingExecutor:
            def __init__(self, *_, **__):
                pass

            def __enter__(self):
                raise RuntimeError("no multiprocessing")

            def __exit__(self, exc_type, exc, tb):
                return False

        ccm = CreateConcatenationMatrix(args)
        monkeypatch.setattr(ccm_module, "ProcessPoolExecutor", FailingExecutor)
        monkeypatch.setattr(
            ccm,
            "_get_taxa_from_alignment",
            lambda path: {path.split("_")[-1]},
        )

        paths = [f"gene_{i}" for i in range(11)]
        taxa = ccm.get_taxa_names(paths)
        assert len(taxa) == 11

    def test_create_missing_seq_str(self, args):
        ccm = CreateConcatenationMatrix(args)
        missing_seq, og_len = ccm.create_missing_seq_str([SeqRecord(Seq("ACGT"), id="A")])
        assert missing_seq == "????"
        assert og_len == 4

    def test_create_missing_seq_str_empty_exits(self, args):
        ccm = CreateConcatenationMatrix(args)
        with pytest.raises(SystemExit) as exc:
            ccm.create_missing_seq_str([])
        assert exc.value.code == 2

    def test_process_taxa_sequences_adds_missing(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="A")]
        concatenated_seqs = defaultdict(list)
        ccm.process_taxa_sequences(records, ["A", "B"], concatenated_seqs, "????")
        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_add_to_partition_info(self, args):
        ccm = CreateConcatenationMatrix(args)
        partition_info, first_len, second_len = ccm.add_to_partition_info([], 10, "AUTO", "g1.fa", 1, 0)
        assert partition_info == ["AUTO, g1.fa=1-10\n"]
        assert first_len == 11
        assert second_len == 10

    def test_add_to_occupancy_info(self, args):
        ccm = CreateConcatenationMatrix(args)
        occ = ccm.add_to_occupancy_info([], {"A"}, ["A", "B"], "g1.fa")
        assert occ == ["g1.fa\t1\t1\t0.5000\tB\n"]

    def test_fasta_and_text_file_writers(self, tmp_path, args):
        ccm = CreateConcatenationMatrix(args)
        fasta_out = tmp_path / "out.fa"
        ccm.fasta_file_write(str(fasta_out), {"A": ["AA", "CC"], "B": ["GG"]})
        assert fasta_out.read_text() == ">A\nAACC\n>B\nGG\n"

        text_out = tmp_path / "out.txt"
        ccm.write_occupancy_or_partition_file(["x\n", "y\n"], str(text_out))
        assert text_out.read_text() == "x\ny\n"

    def test_process_alignment_file_populates_missing_taxa(self, tmp_path):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("C", "AAGT")])
        _, seq_dict, present_taxa, og_len = CreateConcatenationMatrix._process_alignment_file(
            str(fasta),
            ["A", "B", "C"],
        )
        assert present_taxa == {"A", "C"}
        assert og_len == 4
        assert seq_dict["A"] == "ACGT"
        assert seq_dict["B"] == "????"
        assert seq_dict["C"] == "AAGT"

    def test_create_concatenation_matrix_sequential(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "ACGT"), ("B", "A-GT")])
        _write_fasta(gene2, [("A", "TT"), ("C", "TA")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "nested" / "concat"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        occupancy_text = Path(f"{prefix}.occupancy").read_text()
        partition_text = Path(f"{prefix}.partition").read_text()

        assert ">A\nACGTTT\n" in fasta_text
        assert ">B\nA-GT??\n" in fasta_text
        assert ">C\n????TA\n" in fasta_text
        assert "AUTO, " in partition_text
        assert str(gene1) in occupancy_text
        assert str(gene2) in occupancy_text

    def test_create_concatenation_matrix_parallel_fallback(self, tmp_path, monkeypatch):
        class FailingExecutor:
            def __init__(self, *_, **__):
                pass

            def __enter__(self):
                raise RuntimeError("disabled")

            def __exit__(self, exc_type, exc, tb):
                return False

        gene_files = []
        for idx, entries in enumerate(
            [
                [("A", "AA"), ("B", "CC")],
                [("A", "GG"), ("C", "TT")],
                [("B", "TA"), ("C", "AT")],
            ]
        ):
            gene = tmp_path / f"g{idx}.fa"
            _write_fasta(gene, entries)
            gene_files.append(gene)

        alignment_list = tmp_path / "alignments_parallel.txt"
        alignment_list.write_text("\n".join(str(p) for p in gene_files) + "\n")
        prefix = tmp_path / "concat_parallel"

        monkeypatch.setattr(ccm_module, "ProcessPoolExecutor", FailingExecutor)
        ccm = CreateConcatenationMatrix(
            Namespace(alignment_list=str(alignment_list), prefix=str(prefix), json=False, plot_occupancy=False)
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        assert Path(f"{prefix}.fa").exists()
        assert Path(f"{prefix}.occupancy").exists()
        assert Path(f"{prefix}.partition").exists()

    def test_plot_concatenation_occupancy(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.png"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_create_concatenation_matrix_json_and_plot(self, tmp_path, monkeypatch):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "AC"), ("B", "TT")])
        _write_fasta(gene2, [("A", "GA"), ("C", "TA")])
        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_json"
        captured = {}

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=True,
                plot_occupancy=True,
                plot_output=str(tmp_path / "custom_plot.png"),
            )
        )

        monkeypatch.setattr(ccm, "_plot_concatenation_occupancy", lambda **kwargs: captured.setdefault("plot", kwargs))
        monkeypatch.setattr(ccm_module, "print_json", lambda payload: captured.setdefault("payload", payload))
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        assert captured["plot"]["output_file"] == str(tmp_path / "custom_plot.png")
        assert captured["payload"]["total_alignments"] == 2
        assert captured["payload"]["total_taxa"] == 3
        assert captured["payload"]["concatenated_length"] == 4
        assert "occupancy_plot" in captured["payload"]["output_files"]

    def test_compute_effective_occupancy_all_valid(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["ACGT", "TTGG"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 1.0
        assert excluded == set()

    def test_compute_effective_occupancy_all_gaps(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["----", "NNNN"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 0.0
        assert "A" in excluded

    def test_compute_effective_occupancy_mixed(self, args):
        ccm = CreateConcatenationMatrix(args)
        # 4 informative out of 8 total = 0.5
        seqs = {"A": ["AC--", "??GT"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 0.5
        # 0.5 is not < 0.5, so not excluded
        assert excluded == set()

    def test_threshold_excludes_low_occupancy_taxon(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        # A is fully present in both genes
        # B has only gaps/ambiguous chars
        _write_fasta(gene1, [("A", "ACGT"), ("B", "----")])
        _write_fasta(gene2, [("A", "TTGG"), ("B", "NNNN")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_thresh"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=0.5,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" not in fasta_text

    def test_threshold_zero_disables_filtering(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "ACGT"), ("B", "----")])
        _write_fasta(gene2, [("A", "TTGG"), ("B", "NNNN")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_nofilt"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=0,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" in fasta_text

    def test_threshold_one_excludes_all_but_complete(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        # A: fully present in both → occupancy 1.0
        # B: present in gene1, missing from gene2 → has ???? → occupancy 0.5
        _write_fasta(gene1, [("A", "ACGT"), ("B", "TTGG")])
        _write_fasta(gene2, [("A", "AACC")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_strict"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=1.0,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" not in fasta_text
