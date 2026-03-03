"""Lazy service factory declarations for CLI command handlers."""

import importlib
from typing import Dict


class _LazyServiceFactory:
    def __init__(self, module_path: str, class_name: str):
        self.module_path = module_path
        self.class_name = class_name
        self._klass = None

    def __call__(self, *args, **kwargs):
        if self._klass is None:
            module = importlib.import_module(self.module_path)
            self._klass = getattr(module, self.class_name)
        return self._klass(*args, **kwargs)


# Alignment service loaders
AlignmentLength = _LazyServiceFactory("phykit.services.alignment.alignment_length", "AlignmentLength")
AlignmentLengthNoGaps = _LazyServiceFactory("phykit.services.alignment.alignment_length_no_gaps", "AlignmentLengthNoGaps")
AlignmentEntropy = _LazyServiceFactory("phykit.services.alignment.alignment_entropy", "AlignmentEntropy")
AlignmentRecoding = _LazyServiceFactory("phykit.services.alignment.alignment_recoding", "AlignmentRecoding")
AlignmentOutlierTaxa = _LazyServiceFactory("phykit.services.alignment.alignment_outlier_taxa", "AlignmentOutlierTaxa")
ColumnScore = _LazyServiceFactory("phykit.services.alignment.column_score", "ColumnScore")
CompositionalBiasPerSite = _LazyServiceFactory("phykit.services.alignment.compositional_bias_per_site", "CompositionalBiasPerSite")
CompositionPerTaxon = _LazyServiceFactory("phykit.services.alignment.composition_per_taxon", "CompositionPerTaxon")
CreateConcatenationMatrix = _LazyServiceFactory("phykit.services.alignment.create_concatenation_matrix", "CreateConcatenationMatrix")
DNAThreader = _LazyServiceFactory("phykit.services.alignment.dna_threader", "DNAThreader")
EvolutionaryRatePerSite = _LazyServiceFactory("phykit.services.alignment.evolutionary_rate_per_site", "EvolutionaryRatePerSite")
Faidx = _LazyServiceFactory("phykit.services.alignment.faidx", "Faidx")
GCContent = _LazyServiceFactory("phykit.services.alignment.gc_content", "GCContent")
MaskAlignment = _LazyServiceFactory("phykit.services.alignment.mask_alignment", "MaskAlignment")
PlotAlignmentQC = _LazyServiceFactory("phykit.services.alignment.plot_alignment_qc", "PlotAlignmentQC")
OccupancyPerTaxon = _LazyServiceFactory("phykit.services.alignment.occupancy_per_taxon", "OccupancyPerTaxon")
PairwiseIdentity = _LazyServiceFactory("phykit.services.alignment.pairwise_identity", "PairwiseIdentity")
ParsimonyInformative = _LazyServiceFactory("phykit.services.alignment.parsimony_informative_sites", "ParsimonyInformative")
RelativeCompositionVariability = _LazyServiceFactory("phykit.services.alignment.rcv", "RelativeCompositionVariability")
RelativeCompositionVariabilityTaxon = _LazyServiceFactory("phykit.services.alignment.rcvt", "RelativeCompositionVariabilityTaxon")
RenameFastaEntries = _LazyServiceFactory("phykit.services.alignment.rename_fasta_entries", "RenameFastaEntries")
SumOfPairsScore = _LazyServiceFactory("phykit.services.alignment.sum_of_pairs_score", "SumOfPairsScore")
VariableSites = _LazyServiceFactory("phykit.services.alignment.variable_sites", "VariableSites")

AncestralReconstruction = _LazyServiceFactory("phykit.services.tree.ancestral_reconstruction", "AncestralReconstruction")
ConcordanceAsr = _LazyServiceFactory("phykit.services.tree.concordance_asr", "ConcordanceAsr")

# Tree service loaders
BipartitionSupportStats = _LazyServiceFactory("phykit.services.tree.bipartition_support_stats", "BipartitionSupportStats")
BranchLengthMultiplier = _LazyServiceFactory("phykit.services.tree.branch_length_multiplier", "BranchLengthMultiplier")
CollapseBranches = _LazyServiceFactory("phykit.services.tree.collapse_branches", "CollapseBranches")
CovaryingEvolutionaryRates = _LazyServiceFactory("phykit.services.tree.covarying_evolutionary_rates", "CovaryingEvolutionaryRates")
ConsensusNetwork = _LazyServiceFactory("phykit.services.tree.consensus_network", "ConsensusNetwork")
ConsensusTree = _LazyServiceFactory("phykit.services.tree.consensus_tree", "ConsensusTree")
DVMC = _LazyServiceFactory("phykit.services.tree.dvmc", "DVMC")
EvolutionaryRate = _LazyServiceFactory("phykit.services.tree.evolutionary_rate", "EvolutionaryRate")
HiddenParalogyCheck = _LazyServiceFactory("phykit.services.tree.hidden_paralogy_check", "HiddenParalogyCheck")
InternalBranchStats = _LazyServiceFactory("phykit.services.tree.internal_branch_stats", "InternalBranchStats")
InternodeLabeler = _LazyServiceFactory("phykit.services.tree.internode_labeler", "InternodeLabeler")
LastCommonAncestorSubtree = _LazyServiceFactory("phykit.services.tree.last_common_ancestor_subtree", "LastCommonAncestorSubtree")
LBScore = _LazyServiceFactory("phykit.services.tree.lb_score", "LBScore")
MonophylyCheck = _LazyServiceFactory("phykit.services.tree.monophyly_check", "MonophylyCheck")
NearestNeighborInterchange = _LazyServiceFactory("phykit.services.tree.nearest_neighbor_interchange", "NearestNeighborInterchange")
PatristicDistances = _LazyServiceFactory("phykit.services.tree.patristic_distances", "PatristicDistances")
PhylogeneticSignal = _LazyServiceFactory("phykit.services.tree.phylogenetic_signal", "PhylogeneticSignal")
PhylogeneticOrdination = _LazyServiceFactory("phykit.services.tree.phylogenetic_ordination", "PhylogeneticOrdination")
Phylomorphospace = _LazyServiceFactory("phykit.services.tree.phylomorphospace", "Phylomorphospace")
PhylogeneticRegression = _LazyServiceFactory("phykit.services.tree.phylogenetic_regression", "PhylogeneticRegression")
PhylogeneticGLM = _LazyServiceFactory("phykit.services.tree.phylogenetic_glm", "PhylogeneticGLM")
StochasticCharacterMap = _LazyServiceFactory("phykit.services.tree.stochastic_character_map", "StochasticCharacterMap")
ContMap = _LazyServiceFactory("phykit.services.tree.cont_map", "ContMap")
DensityMap = _LazyServiceFactory("phykit.services.tree.density_map", "DensityMap")
Phenogram = _LazyServiceFactory("phykit.services.tree.phenogram", "Phenogram")
Cophylo = _LazyServiceFactory("phykit.services.tree.cophylo", "Cophylo")
RateHeterogeneity = _LazyServiceFactory("phykit.services.tree.rate_heterogeneity", "RateHeterogeneity")
FitContinuous = _LazyServiceFactory("phykit.services.tree.fit_continuous", "FitContinuous")
OUwie = _LazyServiceFactory("phykit.services.tree.ouwie", "OUwie")
OUShiftDetection = _LazyServiceFactory("phykit.services.tree.ou_shift_detection", "OUShiftDetection")
QuartetNetwork = _LazyServiceFactory("phykit.services.tree.quartet_network", "QuartetNetwork")
NetworkSignal = _LazyServiceFactory("phykit.services.tree.network_signal", "NetworkSignal")
LTT = _LazyServiceFactory("phykit.services.tree.ltt", "LTT")
RelativeRateTest = _LazyServiceFactory("phykit.services.tree.relative_rate_test", "RelativeRateTest")
ThresholdModel = _LazyServiceFactory("phykit.services.tree.threshold_model", "ThresholdModel")
PolytomyTest = _LazyServiceFactory("phykit.services.tree.polytomy_test", "PolytomyTest")
PrintTree = _LazyServiceFactory("phykit.services.tree.print_tree", "PrintTree")
PruneTree = _LazyServiceFactory("phykit.services.tree.prune_tree", "PruneTree")
RenameTreeTips = _LazyServiceFactory("phykit.services.tree.rename_tree_tips", "RenameTreeTips")
RobinsonFouldsDistance = _LazyServiceFactory("phykit.services.tree.rf_distance", "RobinsonFouldsDistance")
RootTree = _LazyServiceFactory("phykit.services.tree.root_tree", "RootTree")
Saturation = _LazyServiceFactory("phykit.services.tree.saturation", "Saturation")
SpuriousSequence = _LazyServiceFactory("phykit.services.tree.spurious_sequence", "SpuriousSequence")
TerminalBranchStats = _LazyServiceFactory("phykit.services.tree.terminal_branch_stats", "TerminalBranchStats")
TipLabels = _LazyServiceFactory("phykit.services.tree.tip_labels", "TipLabels")
TipToTipDistance = _LazyServiceFactory("phykit.services.tree.tip_to_tip_distance", "TipToTipDistance")
TipToTipNodeDistance = _LazyServiceFactory("phykit.services.tree.tip_to_tip_node_distance", "TipToTipNodeDistance")
TotalTreeLength = _LazyServiceFactory("phykit.services.tree.total_tree_length", "TotalTreeLength")
Treeness = _LazyServiceFactory("phykit.services.tree.treeness", "Treeness")
TreenessOverRCV = _LazyServiceFactory("phykit.services.tree.treeness_over_rcv", "TreenessOverRCV")
EvoTempoMap = _LazyServiceFactory("phykit.services.tree.evo_tempo_map", "EvoTempoMap")
DiscordanceAsymmetry = _LazyServiceFactory("phykit.services.tree.discordance_asymmetry", "DiscordanceAsymmetry")
SpectralDiscordance = _LazyServiceFactory("phykit.services.tree.spectral_discordance", "SpectralDiscordance")

SERVICE_FACTORIES: Dict[str, _LazyServiceFactory] = {
    name: value
    for name, value in locals().items()
    if isinstance(value, _LazyServiceFactory)
}
