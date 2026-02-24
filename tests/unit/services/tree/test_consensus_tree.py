from argparse import Namespace

import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.consensus_tree import ConsensusTree


def _write(path, content):
    path.write_text(content)
    return str(path)


class TestConsensusTree:
    def test_process_args_defaults_json_false(self):
        args = Namespace(trees="trees.txt", method="majority", missing_taxa="error")
        parsed = ConsensusTree(args).process_args(args)
        assert parsed["json_output"] is False
        assert parsed["method"] == "majority"

    def test_run_json_with_newick_lines(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,D));\n")

        captured = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="error",
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.consensus_tree.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )

        svc.run()

        payload = captured["payload"]
        assert payload["method"] == "majority"
        assert payload["taxa_in_consensus"] == 4
        assert payload["pruned_to_shared_taxa"] is False
        assert "A" in payload["tree_newick"] and "D" in payload["tree_newick"]

    def test_missing_taxa_error_mode_raises(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="error",
                json=False,
            )
        )

        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_taxa_shared_mode_prunes(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        captured = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="shared",
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.consensus_tree.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )

        svc.run()

        payload = captured["payload"]
        assert payload["pruned_to_shared_taxa"] is True
        assert payload["taxa_in_consensus"] == 3
        assert "D" not in payload["tree_newick"]
        assert "E" not in payload["tree_newick"]

    def test_path_list_input(self, tmp_path, monkeypatch):
        t1 = tmp_path / "t1.tre"
        t2 = tmp_path / "t2.tre"
        list_file = tmp_path / "trees.txt"
        _write(t1, "((A,B),(C,D));\n")
        _write(t2, "((A,B),(C,D));\n")
        _write(list_file, f"{t1.name}\n{t2.name}\n")

        printed = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(list_file),
                method="strict",
                missing_taxa="error",
                json=False,
            )
        )
        monkeypatch.setattr("builtins.print", lambda x: printed.setdefault("v", x))

        svc.run()

        assert "A" in printed["v"] and "D" in printed["v"]
