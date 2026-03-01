#!/usr/bin/env Rscript
# Validate PhyKIT consensus_network split frequencies against phangorn::as.splits
#
# Usage: Rscript validate_consensus_network.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

# Helper: extract non-trivial splits as canonicalized strings with frequencies
get_nontrivial_splits <- function(trees) {
  sp    <- as.splits(trees)
  labs  <- attr(sp, "labels")
  wts   <- attr(sp, "weights")
  confs <- attr(sp, "confidences")
  n_taxa <- length(labs)

  results <- data.frame(split = character(), count = integer(),
                        frequency = numeric(), stringsAsFactors = FALSE)

  for (i in seq_along(sp)) {
    taxa_in <- labs[sp[[i]]]
    n_in <- length(taxa_in)
    if (n_in <= 1 || n_in >= (n_taxa - 1)) next

    complement <- setdiff(labs, taxa_in)
    side_a <- sort(taxa_in)
    side_b <- sort(complement)

    if (length(side_a) > length(side_b)) {
      display <- side_b
    } else if (length(side_a) < length(side_b)) {
      display <- side_a
    } else {
      display <- if (paste(side_a, collapse=",") < paste(side_b, collapse=","))
        side_a else side_b
    }

    results <- rbind(results, data.frame(
      split = paste0("{", paste(display, collapse=", "), "}"),
      count = as.integer(wts[i]),
      frequency = round(confs[i], 4),
      stringsAsFactors = FALSE
    ))
  }
  results <- results[order(-results$frequency, results$split), ]
  rownames(results) <- NULL
  results
}

# ---------------------------------------------------------------
# Test 1: Three 4-taxon trees
# ---------------------------------------------------------------
cat("=== Test 1: 2x((A,B),(C,D)) + 1x((A,C),(B,D)) ===\n")
trees1 <- read.tree(text = c(
  "((A,B),(C,D));",
  "((A,B),(C,D));",
  "((A,C),(B,D));"
))
res1 <- get_nontrivial_splits(trees1)
print(res1)
cat("\n")

# ---------------------------------------------------------------
# Test 2: Five identical trees
# ---------------------------------------------------------------
cat("=== Test 2: 5x((A,B),(C,D)) ===\n")
trees2 <- read.tree(text = rep("((A,B),(C,D));", 5))
res2 <- get_nontrivial_splits(trees2)
print(res2)
cat("\n")

# ---------------------------------------------------------------
# Test 3: Six-taxon trees with conflict
# ---------------------------------------------------------------
cat("=== Test 3: 6-taxon trees with conflict ===\n")
trees3 <- read.tree(text = c(
  "((A,B),((C,D),(E,F)));",
  "((A,B),(C,(D,(E,F))));",
  "((A,B),((C,D),(E,F)));",
  "((A,B),(D,(C,(E,F))));",
  "((A,B),((C,D),(E,F)));",
  "((A,B),(C,(D,(E,F))));",
  "((A,B),((C,D),(E,F)));",
  "((A,B),((C,D),(E,F)));",
  "(((B,C),A),(D,(E,F)));",
  "(((A,B,C)),((D),(E,F)));"
))
res3 <- get_nontrivial_splits(trees3)
print(res3)
cat("\n")

# ---------------------------------------------------------------
# Test 4: CSV-style output for automated comparison
# ---------------------------------------------------------------
cat("=== CSV output for automated comparison ===\n")
cat("test,split,count,frequency\n")
for (i in seq_len(nrow(res1))) {
  cat(sprintf("test1,%s,%d,%.4f\n", res1$split[i], res1$count[i], res1$frequency[i]))
}
for (i in seq_len(nrow(res2))) {
  cat(sprintf("test2,%s,%d,%.4f\n", res2$split[i], res2$count[i], res2$frequency[i]))
}
for (i in seq_len(nrow(res3))) {
  cat(sprintf("test3,%s,%d,%.4f\n", res3$split[i], res3$count[i], res3$frequency[i]))
}
