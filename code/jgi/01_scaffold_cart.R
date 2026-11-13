#!/usr/bin/env Rscript
#
# Prepare scaffold lists for import into the JGI IMG/MER Scaffold Cart
#
# Usage:
#   Rscript 01_scaffold_cart.R <scaffolds.tsv> <output_dir> [n_chunks]
#
# Arguments:
#   scaffolds.tsv - tab-separated table with columns `taxon_oid` (IMG genome/metagenome ID)
#                   and `scaffold_id` (e.g., Ga0194039_1011500); one row per selected scaffold
#                   (i.e., top 10 hits by bit score per query and marker, alignment length > 300 bp,
#                   as obtained from the IMG/MER BLAST search with SSU and LSU of the novel taxa)
#   output_dir    - output directory
#   n_chunks      - number of output files (default: 10, i.e. ~100-200 scaffolds per file)
#
# Output:
#   <output_dir>/scaffold_ids_<N>.txt - lists in the IMG/MER import format:
#     "Scaffold ID" header, then one line per scaffold: "<taxon_oid> assembled <scaffold_id>"
#
# Next steps (manual, IMG/MER web interface; requires a JGI account):
#   1. https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=ScaffSearch&page=searchForm
#   2. Upload each file into the Scaffold Cart
#   3. Select all -> "Export Scaffold Data" job; download the resulting *.tar.gz archives
#   4. Continue with 02_extract_rrna.sh
#
# Requirements: R >= 4.0, data.table

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_scaffold_cart.R <scaffolds.tsv> <output_dir> [n_chunks]")
}
out_dir  <- args[2]
n_chunks <- if (length(args) >= 3) as.integer(args[3]) else 10L

scf <- fread(args[1], sep = "\t", colClasses = "character")
if (!all(c("taxon_oid", "scaffold_id") %in% names(scf))) {
  stop("Input table must contain columns `taxon_oid` and `scaffold_id`")
}

scf <- unique(scf[ , .(taxon_oid, scaffold_id) ])
setorder(scf, taxon_oid, scaffold_id)
scf[ , Export := paste0(taxon_oid, " assembled ", scaffold_id) ]

## Split into chunks of roughly equal size
n_chunks <- min(n_chunks, nrow(scf))
chunks   <- split(scf$Export, cut(seq_len(nrow(scf)), breaks = n_chunks, labels = FALSE))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(chunks)) {
  writeLines(c("Scaffold ID", chunks[[i]]),
    con = file.path(out_dir, paste0("scaffold_ids_", i, ".txt")))
}

message("Scaffolds: ", nrow(scf), "; files written: ", length(chunks))
