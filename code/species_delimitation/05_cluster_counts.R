#!/usr/bin/env Rscript
#
# Estimate the number of putative species per lineage by hierarchical clustering
# of ITS distances at a lineage-specific distance threshold
#
# Usage:
#   Rscript 05_cluster_counts.R <distances.rds> <thresholds.tsv> <output.tsv>
#
# Arguments:
#   distances.rds  - named list of `dist` objects (from 02_distances_to_mega.R)
#   thresholds.tsv - tab-separated table with columns `lineage` and `threshold` (distance, e.g. 0.02 = 98% similarity)
#   output.tsv     - output table
#
# The dendrogram is built with complete linkage (main estimate) and single linkage (for comparison),
# and cut at the given threshold. The number of clusters is the estimate of putative species
# observed in the data.
#
# Requirements: R >= 4.0, data.table

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript 05_cluster_counts.R <distances.rds> <thresholds.tsv> <output.tsv>")
}

DD  <- readRDS(args[1])
THR <- fread(args[2], sep = "\t")

missing <- setdiff(THR$lineage, names(DD))
if (length(missing) > 0) {
  message("WARNING: no distance matrix for: ", paste(missing, collapse = ", "))
}
THR <- THR[ lineage %in% names(DD) ]

count_clusters <- function(d, threshold) {
  hc_complete <- hclust(d, method = "complete")
  hc_single   <- hclust(d, method = "single")
  data.table(
    n_seqs              = attr(d, "Size"),
    n_clusters_complete = uniqueN(cutree(hc_complete, h = threshold)),
    n_clusters_single   = uniqueN(cutree(hc_single,   h = threshold)))
}

RES <- THR[ , count_clusters(DD[[lineage]], threshold), by = .(lineage, threshold) ]

fwrite(RES, file = args[3], sep = "\t")
print(RES)
