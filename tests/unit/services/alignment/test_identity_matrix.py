import json
import os
import builtins
import importlib
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.alignment.identity_matrix import IdentityMatrix
from phykit.services.alignment.identity_matrix import _identical_sequence_identity_value
import phykit.services.alignment.identity_matrix as identity_matrix_module


def _write_alignment(path, seqs):
    """Write a dict of {name: seq} as FASTA."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _make_args(alignment, output, **kwargs):
    defaults = dict(
        alignment=str(alignment),
        output=str(output),
        metric="identity",
        tree=None,
        sort="cluster",
        partition=None,
        json=False,
        fig_width=None,
        fig_height=None,
        dpi=300,
        no_title=False,
        title=None,
        legend_position=None,
        ylabel_fontsize=None,
        xlabel_fontsize=None,
        title_fontsize=None,
        axis_fontsize=None,
        colors=None,
        ladderize=False,
        cladogram=False,
        circular=False,
        color_file=None,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


def test_module_import_does_not_import_scipy_clustering(monkeypatch):
    module_name = "phykit.services.alignment.identity_matrix"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.cluster.hierarchy"
            or name.startswith("scipy.cluster.hierarchy.")
            or name == "scipy.spatial.distance"
            or name.startswith("scipy.spatial.distance.")
        ):
            raise AssertionError(
                "identity_matrix module import should not import SciPy clustering"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


def test_module_import_does_not_import_biopython_tree_or_seqio():
    code = """
import sys
import phykit.services.alignment.identity_matrix as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "hashlib" not in sys.modules
assert "phykit.services.tree.base" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert module._LINKAGE is None
assert module._DENDROGRAM is None
assert module._LEAVES_LIST is None
assert module._SQUAREFORM is None
    """
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_proxy_caches_resolved_attributes():
    lazy_np = identity_matrix_module._LazyNumpy()

    first = lazy_np.count_nonzero
    second = lazy_np.count_nonzero

    assert first is second
    assert lazy_np.__dict__["count_nonzero"] is first


def test_identical_sequence_identity_value_ascii_valid_and_invalid():
    assert _identical_sequence_identity_value("-?NX*nxA") == 1.0
    assert _identical_sequence_identity_value("-?NX*nx") == 0.0


def test_identical_sequence_identity_value_unicode_fallback():
    assert _identical_sequence_identity_value("-?NX*nx\u03a9") == 1.0


def test_all_sequences_identical_does_not_slice_taxa_names():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("identity-matrix sequence scan should not slice")
            return super().__getitem__(key)

    taxa_names = NoSliceList(["a", "b", "c"])
    sequences = {"a": "ACGT", "b": "ACGT", "c": "ACGT"}

    assert identity_matrix_module._all_sequences_identical(sequences, taxa_names) is True


def test_repeated_cluster_wrappers_cache_scipy_imports(monkeypatch):
    previous_linkage = identity_matrix_module._LINKAGE
    previous_leaves_list = identity_matrix_module._LEAVES_LIST
    previous_squareform = identity_matrix_module._SQUAREFORM
    identity_matrix_module._LINKAGE = None
    identity_matrix_module._LEAVES_LIST = None
    identity_matrix_module._SQUAREFORM = None
    original_import = builtins.__import__
    cluster_imports = 0
    distance_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal cluster_imports, distance_imports
        if name == "scipy.cluster.hierarchy":
            cluster_imports += 1
        elif name == "scipy.spatial.distance":
            distance_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        distance_matrix = np.array(
            [
                [0.0, 0.2, 0.3],
                [0.2, 0.0, 0.4],
                [0.3, 0.4, 0.0],
            ]
        )
        condensed = identity_matrix_module.squareform(distance_matrix, checks=False)
        cluster = identity_matrix_module.linkage(condensed, method="average")
        identity_matrix_module.leaves_list(cluster)
        first_cluster_imports = cluster_imports
        first_distance_imports = distance_imports

        condensed = identity_matrix_module.squareform(distance_matrix, checks=False)
        cluster = identity_matrix_module.linkage(condensed, method="average")
        identity_matrix_module.leaves_list(cluster)
    finally:
        identity_matrix_module._LINKAGE = previous_linkage
        identity_matrix_module._LEAVES_LIST = previous_leaves_list
        identity_matrix_module._SQUAREFORM = previous_squareform

    assert first_cluster_imports > 0
    assert first_distance_imports > 0
    assert cluster_imports == first_cluster_imports
    assert distance_imports == first_distance_imports


class TestIdentityMatrixUnit:
    def test_parse_alignment_rejects_unequal_sequence_lengths(
        self, tmp_path, capsys, monkeypatch
    ):
        out_path = tmp_path / "out.png"
        args = _make_args("aln.fa", out_path)
        service = IdentityMatrix(args)

        class Record:
            def __init__(self, record_id, seq):
                self.id = record_id
                self.seq = seq

        monkeypatch.setattr(
            service,
            "get_alignment_and_format",
            lambda: (
                [Record("taxon_A", "ACGT"), Record("taxon_B", "ACG")],
                "fasta",
                False,
            ),
        )

        with pytest.raises(SystemExit) as exc_info:
            service._parse_alignment()

        assert exc_info.value.code == 2
        assert capsys.readouterr().err == (
            "Error: sequences have different lengths. Is this an alignment?\n"
        )

    def test_identity_self_is_one(self, tmp_path):
        """Diagonal entries should all be 1.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        for i in range(len(taxa_names)):
            assert matrix[i, i] == 1.0

    def test_identity_symmetric(self, tmp_path):
        """Matrix[i,j] should equal matrix[j,i]."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        n = len(taxa_names)
        for i in range(n):
            for j in range(n):
                assert matrix[i, j] == matrix[j, i]

    def test_identical_sequences(self, tmp_path):
        """Two identical sequences should have identity 1.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        assert matrix[0, 1] == 1.0

    def test_completely_different(self, tmp_path):
        """No matching non-ambiguous positions should give identity 0.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "AAAA",
            "taxon_B": "TTTT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        assert matrix[0, 1] == 0.0

    def test_p_distance_metric(self, tmp_path):
        """p-distance should be 1 - identity."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT",
            "taxon_B": "ACTT",
        })
        args = _make_args(aln_path, out_path, metric="p-distance")
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        identity_matrix = service._compute_identity_matrix(sequences, taxa_names)
        # identity should be 3/4 = 0.75
        assert abs(identity_matrix[0, 1] - 0.75) < 1e-9
        # p-distance = 1 - 0.75 = 0.25
        p_dist_matrix = 1.0 - identity_matrix
        assert abs(p_dist_matrix[0, 1] - 0.25) < 1e-9

    def test_gaps_excluded(self, tmp_path):
        """Gaps and ambiguous chars should be skipped in comparison."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        # Positions: A=A (match), C=C (match), -=G (gap in A, skip), T=T (match)
        # Valid positions: 3, matches: 3, identity = 1.0
        _write_alignment(aln_path, {
            "taxon_A": "AC-T",
            "taxon_B": "ACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        # Position 3 has a gap in taxon_A, so only 3 valid positions, all match
        assert matrix[0, 1] == 1.0

    def test_identity_matrix_fast_path_matches_pairwise_fallback(self, tmp_path):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT-?NXac",
            "taxon_B": "ACCTT?NXat",
            "taxon_C": "TCGTACGTac",
            "taxon_D": "ACGTACGTtt",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()

        fast = service._compute_identity_matrix(sequences, taxa_names)
        legacy = service._compute_identity_matrix_pairwise(sequences, taxa_names)

        np.testing.assert_allclose(fast, legacy)

    def test_identity_matrix_pairwise_fallback_counts_matches_with_count_nonzero(
        self, tmp_path, mocker
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        sequences = {
            "taxon_A": "ACGTΩ-",
            "taxon_B": "ACCTΩ?",
            "taxon_C": "TCGTΩN",
        }
        taxa_names = list(sequences)
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        count_nonzero_spy = mocker.spy(identity_matrix_module.np, "count_nonzero")

        matrix = service._compute_identity_matrix_pairwise(sequences, taxa_names)

        assert matrix[0, 1] == pytest.approx(4 / 5)
        assert matrix[0, 2] == pytest.approx(4 / 5)
        assert any(
            call.args[0].dtype == bool
            for call in count_nonzero_spy.call_args_list
        )

    def test_identity_matrix_clean_ascii_direct_path_skips_validity_mask(
        self, tmp_path, mocker, monkeypatch
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTCGT",
            "taxon_C": "TCGTACGA",
        })
        monkeypatch.setattr(identity_matrix_module, "_NO_INVALID_DIRECT_MIN_LENGTH", 1)
        monkeypatch.setattr(
            identity_matrix_module,
            "_NO_INVALID_DIRECT_SHORT_ALIGNMENT_MIN_TAXA",
            3,
        )
        monkeypatch.setattr(identity_matrix_module, "_NO_INVALID_DIRECT_MAX_TAXA", 10)
        mocker.patch.object(
            identity_matrix_module.np,
            "isin",
            side_effect=AssertionError(
                "clean ASCII identity matrix should skip validity-mask construction"
            ),
        )
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()

        matrix = service._compute_identity_matrix(sequences, taxa_names)

        np.testing.assert_allclose(
            matrix,
            np.array(
                [
                    [1.0, 0.875, 0.75],
                    [0.875, 1.0, 0.625],
                    [0.75, 0.625, 1.0],
                ]
            ),
        )

    def test_identity_matrix_clean_ascii_long_alignment_keeps_small_taxa_generic(
        self, tmp_path, mocker, monkeypatch
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTCGT",
        })
        monkeypatch.setattr(identity_matrix_module, "_NO_INVALID_DIRECT_MIN_LENGTH", 1)
        monkeypatch.setattr(identity_matrix_module, "_NO_INVALID_DIRECT_MIN_TAXA", 3)
        monkeypatch.setattr(
            identity_matrix_module,
            "_NO_INVALID_DIRECT_SHORT_ALIGNMENT_MIN_TAXA",
            4,
        )
        monkeypatch.setattr(identity_matrix_module, "_NO_INVALID_DIRECT_MAX_TAXA", 10)
        isin_spy = mocker.spy(identity_matrix_module.np, "isin")
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()

        matrix = service._compute_identity_matrix(sequences, taxa_names)

        assert isin_spy.called
        np.testing.assert_allclose(matrix, np.array([[1.0, 0.875], [0.875, 1.0]]))

    def test_identity_matrix_all_invalid_pairs_are_zero_off_diagonal(self):
        service = IdentityMatrix.__new__(IdentityMatrix)
        matrix = service._compute_identity_matrix(
            {"taxon_A": "----", "taxon_B": "NNXX"},
            ["taxon_A", "taxon_B"],
        )

        np.testing.assert_allclose(matrix, np.array([[1.0, 0.0], [0.0, 1.0]]))

    def test_identity_matrix_identical_sequences_skip_matrix_setup(self, mocker):
        service = IdentityMatrix.__new__(IdentityMatrix)
        mocker.patch.object(
            identity_matrix_module.np,
            "frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        matrix = service._compute_identity_matrix(
            {
                "taxon_A": "ACGTNN--",
                "taxon_B": "ACGTNN--",
                "taxon_C": "ACGTNN--",
            },
            ["taxon_A", "taxon_B", "taxon_C"],
        )

        np.testing.assert_allclose(np.ones((3, 3), dtype=np.float64), matrix)

    def test_identity_matrix_identical_all_invalid_sequences_keep_zero_off_diagonal(
        self, mocker
    ):
        service = IdentityMatrix.__new__(IdentityMatrix)
        mocker.patch.object(
            identity_matrix_module.np,
            "frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        matrix = service._compute_identity_matrix(
            {
                "taxon_A": "NN--",
                "taxon_B": "NN--",
                "taxon_C": "NN--",
            },
            ["taxon_A", "taxon_B", "taxon_C"],
        )

        np.testing.assert_allclose(
            matrix,
            np.eye(3, dtype=np.float64),
        )

    def test_identity_matrix_handles_empty_sequences(self):
        service = IdentityMatrix.__new__(IdentityMatrix)
        matrix = service._compute_identity_matrix(
            {"taxon_A": "", "taxon_B": ""},
            ["taxon_A", "taxon_B"],
        )

        np.testing.assert_allclose(matrix, np.array([[1.0, 0.0], [0.0, 1.0]]))

    def test_parse_partitions_accepts_raxml_rows_and_skips_non_partitions(
        self, tmp_path
    ):
        partition_file = tmp_path / "test.partition"
        partition_file.write_text(
            "   # ignored\n"
            "\n"
            "AUTO, gene1=1-10\n"
            "DNA,   gene2   =   11   -   20\n"
            "gene3=21-30\n"
            "not a partition\n"
        )
        service = IdentityMatrix.__new__(IdentityMatrix)

        assert service._parse_partitions(str(partition_file)) == [
            ("gene1", 0, 10),
            ("gene2", 10, 20),
            ("gene3", 20, 30),
        ]

    def test_partition_identities_skip_ambiguous_sites(self, tmp_path):
        """Per-partition identity should reuse the same ambiguity rules as the full matrix."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "AC-TAA",
            "taxon_B": "ACGTTT",
            "taxon_C": "NNGTAT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()

        part_names, identities = service._compute_partition_identities(
            sequences,
            taxa_names,
            [("left", 0, 3), ("right", 3, 6)],
        )

        assert part_names == ["left", "right"]
        # Pair order is A-B, A-C, B-C. Ambiguous/gap sites are skipped.
        np.testing.assert_allclose(
            identities,
            np.array([
                [1.0, 1 / 3],
                [0.0, 2 / 3],
                [1.0, 2 / 3],
            ]),
        )

    def test_partition_identity_pairwise_fallback_counts_with_count_nonzero(
        self, tmp_path, mocker
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        sequences = {
            "taxon_A": "ACGTΩ-",
            "taxon_B": "ACCTΩ?",
            "taxon_C": "TCGTΩN",
        }
        taxa_names = list(sequences)
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        count_nonzero_spy = mocker.spy(identity_matrix_module.np, "count_nonzero")

        part_names, identities = service._compute_partition_identities_pairwise(
            sequences,
            taxa_names,
            [("all", 0, 6)],
        )

        assert part_names == ["all"]
        np.testing.assert_allclose(
            identities[:, 0],
            np.array([4 / 5, 4 / 5, 3 / 5]),
        )
        assert any(
            call.args[0].dtype == bool
            for call in count_nonzero_spy.call_args_list
        )

    def test_partition_identities_fast_path_matches_pairwise_fallback(
        self, tmp_path, monkeypatch
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT-?NXac",
            "taxon_B": "ACCTT?NXat",
            "taxon_C": "TCGTACGTac",
            "taxon_D": "ACGTACGTtt",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()
        partitions = [("p1", 0, 4), ("p2", 4, 8), ("p3", 8, 10)]

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError("partition identities should use condensed vectors")

        monkeypatch.setattr(
            identity_matrix_module.np,
            "triu_indices",
            fail_triu_indices,
        )

        fast_names, fast = service._compute_partition_identities(
            sequences, taxa_names, partitions
        )
        legacy_names, legacy = service._compute_partition_identities_pairwise(
            sequences, taxa_names, partitions
        )

        assert fast_names == legacy_names
        np.testing.assert_allclose(fast, legacy)

    def test_partition_identities_identical_sequences_skip_matrix_setup(self, mocker):
        service = IdentityMatrix.__new__(IdentityMatrix)
        mocker.patch.object(
            identity_matrix_module.np,
            "frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        part_names, identities = service._compute_partition_identities(
            {
                "taxon_A": "ACGTNN--",
                "taxon_B": "ACGTNN--",
                "taxon_C": "ACGTNN--",
            },
            ["taxon_A", "taxon_B", "taxon_C"],
            [("valid", 0, 4), ("invalid", 4, 8)],
        )

        assert part_names == ["valid", "invalid"]
        np.testing.assert_allclose(
            identities,
            np.array(
                [
                    [1.0, 0.0],
                    [1.0, 0.0],
                    [1.0, 0.0],
                ],
                dtype=np.float64,
            ),
        )

    def test_partition_identity_strip_matches_pairwise_means(self, monkeypatch):
        n_taxa = 4
        n_parts = 3
        order = [2, 0, 3, 1]
        pairs = [(i, j) for i in range(n_taxa) for j in range(i + 1, n_taxa)]
        part_identities = np.array([
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
            [0.7, 0.8, 0.9],
            [1.0, 0.9, 0.8],
            [0.6, 0.4, 0.2],
            [0.3, 0.2, 0.1],
        ])

        expected = np.zeros((n_taxa, n_parts), dtype=np.float64)
        for row_pos, orig_i in enumerate(order):
            row_pairs = [
                pair_idx
                for pair_idx, pair in enumerate(pairs)
                if orig_i in pair
            ]
            expected[row_pos] = np.mean(part_identities[row_pairs], axis=0)

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError("partition strip should scan condensed rows directly")

        def fail_add_at(*_args, **_kwargs):
            raise AssertionError("partition strip should avoid scatter-style add.at")

        class FailAdd:
            at = staticmethod(fail_add_at)

        monkeypatch.setattr(
            identity_matrix_module.np,
            "triu_indices",
            fail_triu_indices,
        )
        monkeypatch.setattr(identity_matrix_module.np, "add", FailAdd)

        observed = IdentityMatrix._partition_identity_strip(
            part_identities, order, n_taxa, n_parts
        )

        np.testing.assert_allclose(observed, expected)

    def test_summarize_identity_matrix_uses_condensed_distances(self, monkeypatch):
        identity_matrix = np.array(
            [
                [1.0, 0.8, 0.95, 0.4],
                [0.8, 1.0, 0.7, 0.2],
                [0.95, 0.7, 1.0, 0.6],
                [0.4, 0.2, 0.6, 1.0],
            ],
            dtype=float,
        )

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError("identity summary should avoid triangular index arrays")

        def fail_generic_reduction(*_args, **_kwargs):
            raise AssertionError("identity summary should use condensed-vector methods")

        monkeypatch.setattr(
            identity_matrix_module.np,
            "triu_indices",
            fail_triu_indices,
        )
        monkeypatch.setattr(identity_matrix_module.np, "mean", fail_generic_reduction)
        monkeypatch.setattr(identity_matrix_module.np, "min", fail_generic_reduction)
        monkeypatch.setattr(identity_matrix_module.np, "max", fail_generic_reduction)
        monkeypatch.setattr(identity_matrix_module.np, "argmin", fail_generic_reduction)
        monkeypatch.setattr(identity_matrix_module.np, "argmax", fail_generic_reduction)

        observed = IdentityMatrix._summarize_identity_matrix(
            identity_matrix,
            ["taxon_A", "taxon_B", "taxon_C", "taxon_D"],
        )

        assert observed == (
            pytest.approx(0.6083333333333334),
            0.2,
            ["taxon_B", "taxon_D"],
            0.95,
            ["taxon_A", "taxon_C"],
        )

    def test_sort_alpha(self, tmp_path):
        """Alphabetical ordering should sort taxa by name."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "Zebra": "ACGT",
            "Apple": "ACGT",
            "Mango": "ACGT",
        })
        args = _make_args(aln_path, out_path, sort="alpha")
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        order, ordered_labels = service._determine_order(matrix, taxa_names, len(taxa_names))
        assert ordered_labels == ["Apple", "Mango", "Zebra"]

    def test_sort_tree(self, tmp_path):
        """Tree-guided ordering should follow tree tip order."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        tree_path = tmp_path / "tree.nwk"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
            "C": "TTTT",
        })
        with open(tree_path, "w") as f:
            f.write("((C,A),B);")
        args = _make_args(aln_path, out_path, sort="tree", tree=str(tree_path))
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        order, ordered_labels = service._determine_order(matrix, taxa_names, len(taxa_names))
        assert ordered_labels == ["C", "A", "B"]

    def test_cluster_plot_reuses_existing_linkage(self, tmp_path, monkeypatch):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
            "C": "TTTT",
            "D": "TCGA",
        })
        args = _make_args(aln_path, out_path, sort="cluster")
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        order, ordered_labels, cluster_linkage = service._determine_order(
            matrix,
            taxa_names,
            len(taxa_names),
            return_linkage=True,
        )

        def fail_linkage(*_args, **_kwargs):
            raise AssertionError("linkage should be reused")

        monkeypatch.setattr(
            "phykit.services.alignment.identity_matrix.linkage",
            fail_linkage,
        )

        service._plot_heatmap(
            matrix,
            taxa_names,
            order,
            ordered_labels,
            len(taxa_names),
            cluster_linkage=cluster_linkage,
        )

        assert out_path.exists()

    def test_sort_tree_uses_direct_terminal_name_traversal(
        self, tmp_path, monkeypatch
    ):
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        tree_path = tmp_path / "tree.nwk"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
            "C": "TTTT",
        })
        tree_path.write_text("((C,A),B);\n")
        args = _make_args(aln_path, out_path, sort="tree", tree=str(tree_path))
        service = IdentityMatrix(args)
        sequences, taxa_names, _ = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        _, ordered_labels = service._determine_order(
            matrix,
            taxa_names,
            len(taxa_names),
        )

        assert ordered_labels == ["C", "A", "B"]

    def test_sort_tree_requires_tree(self, tmp_path):
        """--sort tree without --tree should raise SystemExit."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
        })
        args = _make_args(aln_path, out_path, sort="tree", tree=None)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        with pytest.raises(SystemExit):
            service._determine_order(matrix, taxa_names, len(taxa_names))

    def test_print_text_output_batches_summary(self, mocker):
        service = IdentityMatrix.__new__(IdentityMatrix)
        service.alignment_file_path = "alignment.fa"
        service.metric = "identity"
        service.sort_method = "cluster"
        service.output_path = "identity.png"
        printed = mocker.patch("builtins.print")

        service._print_text_output(
            12,
            900,
            0.87654321,
            0.123456,
            ["taxon_A", "taxon_B"],
            1.0,
            ["taxon_C", "taxon_D"],
        )

        printed.assert_called_once_with(
            "Sequence Identity Matrix\n"
            "Alignment: alignment.fa\n"
            "Taxa: 12\n"
            "Alignment length: 900\n"
            "Metric: identity\n"
            "Sort: cluster\n"
            "Mean pairwise identity: 0.8765\n"
            "Min: 0.1235 (taxon_A vs taxon_B)\n"
            "Max: 1.0 (taxon_C vs taxon_D)\n"
            "Output: identity.png"
        )

    def test_creates_png(self, tmp_path):
        """The run method should create a non-empty output file."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        service.run()
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_json_output(self, tmp_path, capsys):
        """JSON output should contain correct structure and values."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT",
            "taxon_B": "ACGT",
            "taxon_C": "TTTT",
        })
        args = _make_args(aln_path, out_path, json=True)
        service = IdentityMatrix(args)
        service.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["n_taxa"] == 3
        assert payload["alignment_length"] == 4
        assert payload["metric"] == "identity"
        assert "mean_identity" in payload
        assert "min_identity" in payload
        assert "max_identity" in payload
        assert "taxa_order" in payload
        assert "output_file" in payload
        assert len(payload["taxa_order"]) == 3
        # taxon_A and taxon_B are identical, so max should be 1.0
        assert payload["max_identity"]["value"] == 1.0
