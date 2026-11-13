#!/usr/bin/env Rscript
#
# Convert VSEARCH pairwise identities to distance matrices
#
# Usage:
#   Rscript 02_distances_to_mega.R <dist_dir> <output_dir>
#
# Arguments:
#   dist_dir   - directory with <lineage>.blast6.gz files (from 01_pairwise_distances.sh)
#   output_dir - output directory
#
# Output:
#   <output_dir>/<lineage>.csv  - lower-triangular distance matrix in MEGA CSV format (input for ASAPy)
#   <output_dir>/distances.rds  - named list of `dist` objects (input for 05_cluster_counts.R)
#
# Distance = (100 - percent identity) / 100.
#
# Requirements: R >= 4.0, data.table

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 02_distances_to_mega.R <dist_dir> <output_dir>")
}
dist_dir <- args[1]
out_dir  <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## Replace characters that break ASAPy CSV parsing (commas, line breaks, ()<>[])
sanitize_labels <- function(labels) {
  labels <- trimws(labels)
  labels <- gsub("[,\r\n]", "_", labels, perl = TRUE)
  chartr("()<>[]", "______", labels)
}

## Validate sequence labels
validate_labels <- function(labels, size) {
  labels <- sanitize_labels(as.character(labels))
  if (length(labels) != size)    { stop("Number of labels does not match the matrix size") }
  if (anyNA(labels) || any(!nzchar(labels))) { stop("Empty or missing labels") }
  if (anyDuplicated(labels) > 0) { stop("Labels must be unique after sanitizing") }
  return(labels)
}

## Write a `dist` object in MEGA CSV format:
##   lower-triangular, comma-separated, first row and first column contain labels,
##   the top-left cell is empty
dist_to_megax <- function(x, file, digits = 4) {

  size <- attr(x, "Size")
  if (size < 3L)                { stop("ASAPy requires at least 3 sequences") }
  if (any(!is.finite(unclass(x))) || any(unclass(x) < 0)) {
    stop("Distances must be finite and non-negative")
  }

  labels <- validate_labels(attr(x, "Labels"), size)
  full   <- as.matrix(x)
  output <- matrix("", nrow = size + 1L, ncol = size + 1L)

  output[1L, ] <- c("", labels)
  for (i in seq_len(size)) {
    output[i + 1L, 1L] <- labels[i]
    if (i > 1L) {
      output[i + 1L, 2L:i] <- formatC(full[i, seq_len(i - 1L)], format = "f", digits = digits)
    }
  }

  write.table(output, file = file,
    sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n")
}

## Read VSEARCH blast6 output and convert to a `dist` object
vsearch_to_dist <- function(x) {

  bc <- fread(cmd = paste("gzip -dc", shQuote(x)), header = FALSE, sep = "\t",
    col.names = c("Seq1", "Seq2", "Identity", "Length", "Mismatches", "Gaps",
                  "QueryStart", "QueryEnd", "TargetStart", "TargetEnd", "Evalue", "BitScore"))

  bc <- bc[ , .(Seq1, Seq2, Dist = (100 - Identity) / 100) ]

  ## Long table -> square matrix -> dist (lower triangle)
  dd <- with(bc, xtabs(Dist ~ Seq1 + Seq2))
  if (nrow(dd) != ncol(dd)) { message("WARNING: non-square distance matrix") }
  ids <- intersect(rownames(dd), colnames(dd))
  dd  <- as.dist(dd[ids, ids])

  ## Missing pairs end up as zero distances (all input sequences are unique, so a true zero is not possible).
  ## Zero distances make ASAPy crash -> replace them with the maximum distance plus small random noise
  if (any(dd == 0)) {
    nn <- sum(dd == 0)
    message("WARNING: zero distances found (n = ", nn, "), replacing them")
    mx <- max(dd)
    mn <- min(dd[dd > 0])
    set.seed(nn)
    mr <- mx + round(abs(rnorm(n = nn, mean = mn, sd = mn * 2)), 4)
    mr[mr > 1] <- 1
    dd[dd == 0] <- mr
  }

  return(dd)
}


fls <- list.files(path = dist_dir, pattern = "\\.blast6\\.gz$", full.names = TRUE)
if (length(fls) == 0) { stop("No *.blast6.gz files found in ", dist_dir) }

DD <- lapply(fls, function(x) {
  taxon <- sub("\\.blast6\\.gz$", "", basename(x))
  message("..", taxon)
  dd <- vsearch_to_dist(x)
  dist_to_megax(dd, file = file.path(out_dir, paste0(taxon, ".csv")), digits = 4)
  dd
})
names(DD) <- sub("\\.blast6\\.gz$", "", basename(fls))

saveRDS(DD, file = file.path(out_dir, "distances.rds"))
message("Done: ", out_dir)
