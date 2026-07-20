#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: validate_covarying_rates.R <gene_tree_0> <gene_tree_1> <reference_tree>")
}

trees <- lapply(args, read.tree)
shared <- Reduce(intersect, lapply(trees, function(tree) tree$tip.label))
trees <- lapply(trees, keep.tip, tip = shared)

branch_lengths_by_clade <- function(tree) {
  descendants <- vector("list", max(tree$edge))
  for (tip in seq_along(tree$tip.label)) {
    descendants[[tip]] <- tree$tip.label[[tip]]
  }
  for (edge_index in rev(seq_len(nrow(tree$edge)))) {
    parent <- tree$edge[edge_index, 1]
    child <- tree$edge[edge_index, 2]
    descendants[[parent]] <- union(descendants[[parent]], descendants[[child]])
  }

  keys <- vapply(
    tree$edge[, 2],
    function(node) paste(sort(descendants[[node]]), collapse = ";"),
    character(1)
  )
  setNames(tree$edge.length, keys)
}

lengths <- lapply(trees, branch_lengths_by_clade)
reference_keys <- names(lengths[[3]])
rates_zero <- unname(lengths[[1]][reference_keys] / lengths[[3]][reference_keys])
rates_one <- unname(lengths[[2]][reference_keys] / lengths[[3]][reference_keys])

keep <- is.finite(rates_zero) & is.finite(rates_one) &
  abs(rates_zero) <= 5 & abs(rates_one) <= 5
result <- cor.test(rates_zero[keep], rates_one[keep], method = "pearson")

cat(sprintf("%.15f\t%.15f\t%d\n", result$estimate, result$p.value, sum(keep)))
