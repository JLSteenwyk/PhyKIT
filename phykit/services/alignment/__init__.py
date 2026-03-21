import importlib


_EXPORTS = {
    "AlignmentLength": "alignment_length",
    "AlignmentSubsample": "alignment_subsample",
    "AlignmentLengthNoGaps": "alignment_length_no_gaps",
    "AlignmentEntropy": "alignment_entropy",
    "AlignmentRecoding": "alignment_recoding",
    "ColumnScore": "column_score",
    "CompositionalBiasPerSite": "compositional_bias_per_site",
    "CompositionPerTaxon": "composition_per_taxon",
    "CreateConcatenationMatrix": "create_concatenation_matrix",
    "DNAThreader": "dna_threader",
    "Dstatistic": "dstatistic",
    "Dfoil": "dfoil",
    "EvolutionaryRatePerSite": "evolutionary_rate_per_site",
    "Faidx": "faidx",
    "GCContent": "gc_content",
    "MaskAlignment": "mask_alignment",
    "PairwiseIdentity": "pairwise_identity",
    "IdentityMatrix": "identity_matrix",
    "ParsimonyInformative": "parsimony_informative_sites",
    "OccupancyPerTaxon": "occupancy_per_taxon",
    "RelativeCompositionVariability": "rcv",
    "RelativeCompositionVariabilityTaxon": "rcvt",
    "RenameFastaEntries": "rename_fasta_entries",
    "SumOfPairsScore": "sum_of_pairs_score",
    "PhyloGwas": "phylo_gwas",
    "VariableSites": "variable_sites",
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
