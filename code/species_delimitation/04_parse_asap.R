#!/usr/bin/env Rscript
#
# Summarize ASAP results across lineages
#
# Usage:
#   Rscript 04_parse_asap.R <asap_dir> <output_prefix>
#
# Arguments:
#   asap_dir      - directory with one ASAPy output subdirectory per lineage (from 03_run_asap.sh)
#   output_prefix - prefix for output files
#
# Output:
#   <output_prefix>_best_partition.tsv     - lineage, threshold distance and number of species of the top-ranked partition
#   <output_prefix>_scores.tsv             - all reported partitions with their ASAP scores
#   <output_prefix>_partitions_long.tsv    - sequence-to-species assignment in the top-ranked partition
#   <output_prefix>_partitions_wide.tsv    - members (accession IDs) of each putative species
#
# Requirements: R >= 4.0, data.table; sed

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 04_parse_asap.R <asap_dir> <output_prefix>")
}
asap_dir   <- args[1]
out_prefix <- args[2]

drr <- list.dirs(path = asap_dir, full.names = TRUE, recursive = FALSE)

parse_asap <- function(x) {

  ll <- basename(x)
  message("..", ll)

  ## Partition scores (whitespace-formatted text table -> tab-separated)
  sc <- fread(
    cmd = paste0("sed '2d ; /^$/d ; s/^*// ; s/spec w/spec_w/ ; s/^[ \t]*// ; s/[ \t]*$// ; s/ \t/\t/ ; s/#//g ; s/ \\+/\t/g' ",
                 shQuote(file.path(x, "scores.tab"))),
    sep = "\t")

  ## Top-ranked partition (long format)
  pp <- list.files(path = x, pattern = "\\.Partition_1\\.csv$", full.names = TRUE)
  pp <- fread(cmd = paste0("sed 's/ , /\t/' ", shQuote(pp)),
    header = FALSE, col.names = c("SeqName", "Partition"), sep = "\t")
  pp[ , AccID := tstrsplit(SeqName, split = ";", keep = 1) ]

  ## Wide format: one row per putative species
  ww <- pp[ , .(Members = paste(sort(unique(AccID)), collapse = ", ")), by = "Partition"]

  sc[ , Lineage := ll ]; setcolorder(sc, "Lineage")
  pp[ , Lineage := ll ]; setcolorder(pp, "Lineage")
  ww[ , Lineage := ll ]; setcolorder(ww, "Lineage")

  list(
    scores     = sc,
    best_score = sc[ 1, .(Lineage, distance, species) ],
    partition  = pp,
    wide       = ww)
}

ASAPS <- lapply(drr, parse_asap)

fwrite(rbindlist(lapply(ASAPS, `[[`, "best_score")), paste0(out_prefix, "_best_partition.tsv"),  sep = "\t")
fwrite(rbindlist(lapply(ASAPS, `[[`, "scores"), fill = TRUE), paste0(out_prefix, "_scores.tsv"), sep = "\t")
fwrite(rbindlist(lapply(ASAPS, `[[`, "partition")),  paste0(out_prefix, "_partitions_long.tsv"), sep = "\t")
fwrite(rbindlist(lapply(ASAPS, `[[`, "wide")),       paste0(out_prefix, "_partitions_wide.tsv"), sep = "\t")
