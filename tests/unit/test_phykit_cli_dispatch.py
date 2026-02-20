import pytest

import phykit.phykit as phykit_module
from phykit.phykit import Phykit
from phykit.errors import PhykitUserError


COMMAND_METHODS = [
    "alignment_length",
    "alignment_length_no_gaps",
    "alignment_entropy",
    "alignment_recoding",
    "alignment_outlier_taxa",
    "column_score",
    "compositional_bias_per_site",
    "composition_per_taxon",
    "evolutionary_rate_per_site",
    "faidx",
    "gc_content",
    "mask_alignment",
    "plot_alignment_qc",
    "occupancy_per_taxon",
    "pairwise_identity",
    "parsimony_informative_sites",
    "rcv",
    "rcvt",
    "rename_fasta_entries",
    "sum_of_pairs_score",
    "variable_sites",
    "bipartition_support_stats",
    "branch_length_multiplier",
    "collapse_branches",
    "covarying_evolutionary_rates",
    "dvmc",
    "evolutionary_rate",
    "hidden_paralogy_check",
    "internal_branch_stats",
    "internode_labeler",
    "last_common_ancestor_subtree",
    "lb_score",
    "monophyly_check",
    "nearest_neighbor_interchange",
    "patristic_distances",
    "polytomy_test",
    "print_tree",
    "prune_tree",
    "rename_tree_tips",
    "rf_distance",
    "root_tree",
    "spurious_sequence",
    "terminal_branch_stats",
    "tip_labels",
    "tip_to_tip_distance",
    "tip_to_tip_node_distance",
    "total_tree_length",
    "treeness",
    "saturation",
    "treeness_over_rcv",
    "create_concatenation_matrix",
    "thread_dna",
]

COMMAND_FACTORY_METHODS = [
    ("alignment_length", "AlignmentLength"),
    ("alignment_length_no_gaps", "AlignmentLengthNoGaps"),
    ("alignment_entropy", "AlignmentEntropy"),
    ("alignment_recoding", "AlignmentRecoding"),
    ("alignment_outlier_taxa", "AlignmentOutlierTaxa"),
    ("column_score", "ColumnScore"),
    ("compositional_bias_per_site", "CompositionalBiasPerSite"),
    ("composition_per_taxon", "CompositionPerTaxon"),
    ("evolutionary_rate_per_site", "EvolutionaryRatePerSite"),
    ("faidx", "Faidx"),
    ("gc_content", "GCContent"),
    ("mask_alignment", "MaskAlignment"),
    ("plot_alignment_qc", "PlotAlignmentQC"),
    ("occupancy_per_taxon", "OccupancyPerTaxon"),
    ("pairwise_identity", "PairwiseIdentity"),
    ("parsimony_informative_sites", "ParsimonyInformative"),
    ("rcv", "RelativeCompositionVariability"),
    ("rcvt", "RelativeCompositionVariabilityTaxon"),
    ("rename_fasta_entries", "RenameFastaEntries"),
    ("sum_of_pairs_score", "SumOfPairsScore"),
    ("variable_sites", "VariableSites"),
    ("bipartition_support_stats", "BipartitionSupportStats"),
    ("branch_length_multiplier", "BranchLengthMultiplier"),
    ("collapse_branches", "CollapseBranches"),
    ("covarying_evolutionary_rates", "CovaryingEvolutionaryRates"),
    ("dvmc", "DVMC"),
    ("evolutionary_rate", "EvolutionaryRate"),
    ("hidden_paralogy_check", "HiddenParalogyCheck"),
    ("internal_branch_stats", "InternalBranchStats"),
    ("internode_labeler", "InternodeLabeler"),
    ("last_common_ancestor_subtree", "LastCommonAncestorSubtree"),
    ("lb_score", "LBScore"),
    ("monophyly_check", "MonophylyCheck"),
    ("nearest_neighbor_interchange", "NearestNeighborInterchange"),
    ("patristic_distances", "PatristicDistances"),
    ("polytomy_test", "PolytomyTest"),
    ("print_tree", "PrintTree"),
    ("prune_tree", "PruneTree"),
    ("rename_tree_tips", "RenameTreeTips"),
    ("rf_distance", "RobinsonFouldsDistance"),
    ("root_tree", "RootTree"),
    ("spurious_sequence", "SpuriousSequence"),
    ("terminal_branch_stats", "TerminalBranchStats"),
    ("tip_labels", "TipLabels"),
    ("tip_to_tip_distance", "TipToTipDistance"),
    ("tip_to_tip_node_distance", "TipToTipNodeDistance"),
    ("total_tree_length", "TotalTreeLength"),
    ("treeness", "Treeness"),
    ("saturation", "Saturation"),
    ("treeness_over_rcv", "TreenessOverRCV"),
    ("create_concatenation_matrix", "CreateConcatenationMatrix"),
    ("thread_dna", "DNAThreader"),
]


class TestPhykitCliDispatch:
    @pytest.mark.parametrize("method_name", COMMAND_METHODS)
    def test_command_parser_help_exits_cleanly(self, method_name):
        method = getattr(Phykit, method_name)
        with pytest.raises(SystemExit) as exc:
            method(["-h"])
        assert exc.value.code == 0

    def test_run_alias_invalid_command(self, capsys):
        instance = object.__new__(Phykit)
        with pytest.raises(SystemExit) as exc:
            instance.run_alias("not_a_real_cmd", [])
        assert exc.value.code == 1
        out, _ = capsys.readouterr()
        assert "Invalid command option" in out

    def test_run_alias_version(self, capsys):
        instance = object.__new__(Phykit)
        instance.help_header = Phykit.help_header
        instance.run_alias("v", [])
        out, _ = capsys.readouterr()
        assert "Version:" in out

    def test_init_dispatches_named_command(self, monkeypatch):
        calls = {}

        def fake_alignment_length(self, argv):
            calls["argv"] = argv

        monkeypatch.setattr(Phykit, "alignment_length", fake_alignment_length)
        monkeypatch.setattr("sys.argv", ["phykit", "alignment_length", "x.fa"])
        Phykit()
        assert calls["argv"] == ["x.fa"]

    def test_init_dispatches_alias(self, monkeypatch):
        calls = {}

        def fake_run_alias(self, command, argv):
            calls["command"] = command
            calls["argv"] = argv

        monkeypatch.setattr(Phykit, "run_alias", fake_run_alias)
        monkeypatch.setattr("sys.argv", ["phykit", "al", "x.fa"])
        Phykit()
        assert calls == {"command": "al", "argv": ["x.fa"]}

    def test_init_nameerror_exits_2(self, monkeypatch):
        def fake_run_alias(self, command, argv):
            raise NameError("boom")

        monkeypatch.setattr(Phykit, "run_alias", fake_run_alias)
        monkeypatch.setattr("sys.argv", ["phykit", "bad_alias"])
        with pytest.raises(SystemExit) as exc:
            Phykit()
        assert exc.value.code == 2

    def test_main_invokes_phykit(self, monkeypatch):
        calls = {"count": 0}

        def fake_phykit():
            calls["count"] += 1

        monkeypatch.setattr(phykit_module, "Phykit", fake_phykit)
        phykit_module.main()
        assert calls["count"] == 1

    @pytest.mark.parametrize(
        "wrapper_name,method_name",
        [
            ("alignment_length", "alignment_length"),
            ("alignment_length_no_gaps", "alignment_length_no_gaps"),
            ("alignment_entropy", "alignment_entropy"),
            ("alignment_outlier_taxa", "alignment_outlier_taxa"),
            ("column_score", "column_score"),
            ("compositional_bias_per_site", "compositional_bias_per_site"),
            ("composition_per_taxon", "composition_per_taxon"),
            ("evolutionary_rate_per_site", "evolutionary_rate_per_site"),
            ("faidx", "faidx"),
            ("gc_content", "gc_content"),
            ("mask_alignment", "mask_alignment"),
            ("plot_alignment_qc", "plot_alignment_qc"),
            ("occupancy_per_taxon", "occupancy_per_taxon"),
            ("pairwise_identity", "pairwise_identity"),
            ("parsimony_informative_sites", "parsimony_informative_sites"),
            ("rcv", "rcv"),
            ("rcvt", "rcvt"),
            ("rename_fasta_entries", "rename_fasta_entries"),
            ("sum_of_pairs_score", "sum_of_pairs_score"),
            ("variable_sites", "variable_sites"),
            ("bipartition_support_stats", "bipartition_support_stats"),
            ("branch_length_multiplier", "branch_length_multiplier"),
            ("collapse_branches", "collapse_branches"),
            ("covarying_evolutionary_rates", "covarying_evolutionary_rates"),
            ("dvmc", "dvmc"),
            ("evolutionary_rate", "evolutionary_rate"),
            ("hidden_paralogy_check", "hidden_paralogy_check"),
            ("internal_branch_stats", "internal_branch_stats"),
            ("internode_labeler", "internode_labeler"),
            ("last_common_ancestor_subtree", "last_common_ancestor_subtree"),
            ("lb_score", "lb_score"),
            ("monophyly_check", "monophyly_check"),
            ("nearest_neighbor_interchange", "nearest_neighbor_interchange"),
            ("patristic_distances", "patristic_distances"),
            ("polytomy_test", "polytomy_test"),
            ("print_tree", "print_tree"),
            ("prune_tree", "prune_tree"),
            ("rename_tree_tips", "rename_tree_tips"),
            ("rf_distance", "rf_distance"),
            ("root_tree", "root_tree"),
            ("spurious_sequence", "spurious_sequence"),
            ("terminal_branch_stats", "terminal_branch_stats"),
            ("tip_labels", "tip_labels"),
            ("tip_to_tip_distance", "tip_to_tip_distance"),
            ("tip_to_tip_node_distance", "tip_to_tip_node_distance"),
            ("total_tree_length", "total_tree_length"),
            ("treeness", "treeness"),
            ("saturation", "saturation"),
            ("treeness_over_rcv", "treeness_over_rcv"),
            ("create_concatenation_matrix", "create_concatenation_matrix"),
            ("thread_dna", "thread_dna"),
        ],
    )
    def test_module_wrapper_forwards_argv(self, monkeypatch, wrapper_name, method_name):
        calls = {}

        def fake_method(argv):
            calls["argv"] = argv

        monkeypatch.setattr(phykit_module.Phykit, method_name, staticmethod(fake_method))
        monkeypatch.setattr(phykit_module.sys, "argv", ["phykit", "A", "B"])
        getattr(phykit_module, wrapper_name)()
        assert calls["argv"] == ["A", "B"]

    @pytest.mark.parametrize("method_name,factory_name", COMMAND_FACTORY_METHODS)
    def test_command_executes_service_factory(self, monkeypatch, method_name, factory_name):
        calls = {"ran": False}

        class FakeRunner:
            def run(self):
                calls["ran"] = True

        class FakeParser:
            def add_argument(self, *args, **kwargs):
                return None

            def parse_args(self, argv):
                return object()

        monkeypatch.setattr(phykit_module, "_new_parser", lambda **kwargs: FakeParser())
        monkeypatch.setattr(phykit_module, factory_name, lambda args: FakeRunner())

        method = getattr(phykit_module.Phykit, method_name)
        method(["ignored"])
        assert calls["ran"] is True

    def test_run_service_renders_phykit_user_error(self, monkeypatch, capsys):
        class FakeParser:
            def parse_args(self, argv):
                return object()

        class FakeRunner:
            def run(self):
                raise PhykitUserError(["line one", "line two"], code=2)

        with pytest.raises(SystemExit) as exc:
            phykit_module._run_service(FakeParser(), [], lambda args: FakeRunner())

        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "line one" in out
        assert "line two" in out
