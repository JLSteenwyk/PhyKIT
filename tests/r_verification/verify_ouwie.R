#!/usr/bin/env Rscript
# verify_ouwie.R
#
# Verify PhyKIT OUwie implementation against R's OUwie package (v2.10).
# Requires: ape, OUwie, phytools
#
# Usage:
#   Rscript tests/r_verification/verify_ouwie.R
#
# The tree is rooted at "dog" to match BioPython's handling of the
# trifurcating Newick. Regime assignments use ML ancestral state
# reconstruction (ace), which matches PhyKIT's Fitch parsimony for
# this dataset (aquatic node = MRCA of sea_lion + seal; all others
# terrestrial).
#
# Key findings:
#   - BM1: PhyKIT matches R OUwie exactly (root.station=FALSE)
#   - BMS: PhyKIT within 0.07 LL units of R (optimizer convergence)
#   - OU1/OUM/OUMV: R converges to alpha=0 (BM boundary); PhyKIT
#     finds genuinely better OU optima via multi-interval search
#   - OUMA: PhyKIT within 0.003 LL units of R
#   - OUMVA: PhyKIT within 0.02 LL units of R

suppressPackageStartupMessages({
  library(ape)
  library(OUwie)
  library(phytools)
})

# Read and root tree
tree <- read.tree("tests/sample_files/tree_simple.tre")
tree <- root(tree, "dog", resolve.root = TRUE)
tree$edge.length[tree$edge.length == 0] <- 1e-6

# Read data
trait_data <- read.delim("tests/sample_files/tree_simple_traits.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "trait"), stringsAsFactors = FALSE)
regime_data <- read.delim("tests/sample_files/tree_simple_regimes.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "regime"), stringsAsFactors = FALSE)

data <- merge(regime_data, trait_data, by = "taxon")
colnames(data) <- c("species", "regime", "trait")

# Assign regime labels to internal nodes via ML ACE
tip_regimes <- setNames(data$regime, data$species)
regime_vec <- tip_regimes[tree$tip.label]
regime_factor <- as.factor(regime_vec)
names(regime_factor) <- tree$tip.label
anc <- ace(regime_factor, tree, type = "discrete", method = "ML")
node_regimes <- apply(anc$lik.anc, 1, function(x)
  colnames(anc$lik.anc)[which.max(x)])
tree$node.label <- node_regimes

cat("Node labels:", tree$node.label, "\n\n")

# Run all models with root.station=FALSE
# (matches PhyKIT's conditional-on-root VCV)
cat("Model       logLik       AIC     AICc    k\n")
cat("--------------------------------------------\n")

for (model in c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")) {
  tryCatch({
    result <- OUwie(tree, data, model = model, root.station = FALSE,
                    algorithm = "invert", quiet = TRUE,
                    check.identify = FALSE)
    cat(sprintf("%-8s %10.6f %9.4f %9.4f %4d\n",
                model, result$loglik, result$AIC, result$AICc,
                result$param.count))
  }, error = function(e) {
    cat(sprintf("%-8s ERROR: %s\n", model, e$message))
  })
}

# Also produce geiger BM/OU for cross-validation
suppressPackageStartupMessages(library(geiger))
trait_vec <- setNames(data$trait, data$species)
bm <- fitContinuous(tree, trait_vec, model = "BM")
cat(sprintf("\ngeiger BM: logLik=%.6f sigma2=%.6f z0=%.6f\n",
            bm$opt$lnL, bm$opt$sigsq, bm$opt$z0))

cat("\nDone.\n")
