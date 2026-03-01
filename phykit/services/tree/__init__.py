import importlib


_EXPORTS = {
    "BipartitionSupportStats": "bipartition_support_stats",
    "BranchLengthMultiplier": "branch_length_multiplier",
    "CollapseBranches": "collapse_branches",
    "CovaryingEvolutionaryRates": "covarying_evolutionary_rates",
    "ConsensusNetwork": "consensus_network",
    "ConsensusTree": "consensus_tree",
    "DVMC": "dvmc",
    "EvolutionaryRate": "evolutionary_rate",
    "HiddenParalogyCheck": "hidden_paralogy_check",
    "InternalBranchStats": "internal_branch_stats",
    "InternodeLabeler": "internode_labeler",
    "LastCommonAncestorSubtree": "last_common_ancestor_subtree",
    "LBScore": "lb_score",
    "MonophylyCheck": "monophyly_check",
    "NearestNeighborInterchange": "nearest_neighbor_interchange",
    "PatristicDistances": "patristic_distances",
    "PolytomyTest": "polytomy_test",
    "PrintTree": "print_tree",
    "PruneTree": "prune_tree",
    "RenameTreeTips": "rename_tree_tips",
    "RobinsonFouldsDistance": "rf_distance",
    "RootTree": "root_tree",
    "Saturation": "saturation",
    "SpuriousSequence": "spurious_sequence",
    "TerminalBranchStats": "terminal_branch_stats",
    "TipLabels": "tip_labels",
    "TipToTipDistance": "tip_to_tip_distance",
    "TipToTipNodeDistance": "tip_to_tip_node_distance",
    "TotalTreeLength": "total_tree_length",
    "Treeness": "treeness",
    "TreenessOverRCV": "treeness_over_rcv",
}

__all__ = list(_EXPORTS.keys())


def __getattr__(name):
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = importlib.import_module(f".{module_name}", __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value
