import importlib


_EXPORTS = {
    "BipartitionSupportStats": "bipartition_support_stats",
    "BranchLengthMultiplier": "branch_length_multiplier",
    "CollapseBranches": "collapse_branches",
    "CovaryingEvolutionaryRates": "covarying_evolutionary_rates",
    "ConsensusNetwork": "consensus_network",
    "NeighborNet": "neighbor_net",
    "ConsensusTree": "consensus_tree",
    "DVMC": "dvmc",
    "DiscordanceAsymmetry": "discordance_asymmetry",
    "EvolutionaryRate": "evolutionary_rate",
    "EvoTempoMap": "evo_tempo_map",
    "FitDiscrete": "fit_discrete",
    "HiddenParalogyCheck": "hidden_paralogy_check",
    "Hybridization": "hybridization",
    "InternalBranchStats": "internal_branch_stats",
    "InternodeLabeler": "internode_labeler",
    "LastCommonAncestorSubtree": "last_common_ancestor_subtree",
    "LBScore": "lb_score",
    "LTT": "ltt",
    "MonophylyCheck": "monophyly_check",
    "NearestNeighborInterchange": "nearest_neighbor_interchange",
    "PatristicDistances": "patristic_distances",
    "QuartetNetwork": "quartet_network",
    "NetworkSignal": "network_signal",
    "RelativeRateTest": "relative_rate_test",
    "PolytomyTest": "polytomy_test",
    "PrintTree": "print_tree",
    "ParsimonyScore": "parsimony_score",
    "CharacterMap": "character_map",
    "PhyloHeatmap": "phylo_heatmap",
    "PruneTree": "prune_tree",
    "QuartetPie": "quartet_pie",
    "RenameTreeTips": "rename_tree_tips",
    "IndependentContrasts": "independent_contrasts",
    "KuhnerFelsensteinDistance": "kf_distance",
    "RobinsonFouldsDistance": "rf_distance",
    "RootTree": "root_tree",
    "Saturation": "saturation",
    "SimmapSummary": "simmap_summary",
    "Spr": "spr",
    "SpuriousSequence": "spurious_sequence",
    "TerminalBranchStats": "terminal_branch_stats",
    "TipLabels": "tip_labels",
    "TipToTipDistance": "tip_to_tip_distance",
    "TipToTipNodeDistance": "tip_to_tip_node_distance",
    "TotalTreeLength": "total_tree_length",
    "TransferAnnotations": "transfer_annotations",
    "Treeness": "treeness",
    "ThresholdModel": "threshold_model",
    "TreenessOverRCV": "treeness_over_rcv",
    "ConcordanceAsr": "concordance_asr",
    "PhyloLogistic": "phylo_logistic",
    "PhyloAnova": "phylo_anova",
    "PhyloPath": "phylo_path",
    "PhyloImpute": "phylo_impute",
    "TraitCorrelation": "trait_correlation",
    "TraitRateMap": "trait_rate_map",
    "TreeSpace": "tree_space",
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
