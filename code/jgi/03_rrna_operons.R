#!/usr/bin/env Rscript
#
# Extract SSU, ITS and LSU regions from metagenome scaffolds based on Infernal hits
#
# Usage:
#   Rscript 03_rrna_operons.R <input_dir> <output_dir>
#
# Arguments:
#   input_dir  - directory with JGI_scaffolds.fasta, cmsearch_SSU.tbl and cmsearch_LSU.tbl (from 02_extract_rrna.sh)
#   output_dir - output directory
#
# Approach:
#   - Only non-truncated SSU (RF01960) and LSU (RF02543) hits are used
#   - Hits on each scaffold and strand are ordered 5'->3'; an SSU hit followed by an LSU hit
#     on the same strand forms a complete operon; the ITS region is the segment between them
#   - Regions on the minus strand are reverse-complemented
#   - Scaffolds without a complete SSU-ITS-LSU operon are not exported
#
# Output (sequence names: <scaffold>___operon_<N>__<region>):
#   <output_dir>/cmsearch__FullOperon.fasta
#   <output_dir>/cmsearch__SSU.fasta
#   <output_dir>/cmsearch__ITS.fasta
#   <output_dir>/cmsearch__LSU.fasta
#   <output_dir>/cmsearch__operons.tsv  - coordinates of operons and regions
#
# Requirements: R >= 4.0, data.table, Biostrings (Bioconductor)

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 03_rrna_operons.R <input_dir> <output_dir>")
}
inp_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## Read cmsearch --tblout (space-delimited; the last column, description, may contain spaces)
read_tblout <- function(file) {
  lns <- readLines(file)
  lns <- lns[ !grepl("^#", lns) & nzchar(lns) ]
  if (length(lns) == 0) { return(data.table()) }
  x <- tstrsplit(lns, split = " +", keep = c(1, 8, 9, 10, 11, 15, 16))
  data.table(
    scaffold = x[[1]],
    seq_from = as.integer(x[[2]]),
    seq_to   = as.integer(x[[3]]),
    strand   = x[[4]],
    trunc    = x[[5]],
    score    = as.numeric(x[[6]]),
    E_value  = as.numeric(x[[7]]))
}

## Find complete SSU -> LSU operons on one scaffold
find_operons <- function(x) {
  # x = hits for a single scaffold (columns: scaffold, start, end, strand, feature)

  res <- list()
  op_i <- 0L

  for (st in unique(x$strand)) {
    g <- x[ strand == st ]

    ## Biological 5' -> 3' order
    if (st == "-") { setorder(g, -end, -start) } else { setorder(g, start, end) }

    pending_ssu <- NULL
    for (i in seq_len(nrow(g))) {
      if (g$feature[i] == "SSU") {
        pending_ssu <- g[i]          # a new SSU before an LSU replaces the previous (partial) one
      } else if (!is.null(pending_ssu)) {
        lsu <- g[i]
        op_i <- op_i + 1L
        if (st == "-") {
          its_start <- lsu$end + 1L
          its_end   <- pending_ssu$start - 1L
        } else {
          its_start <- pending_ssu$end + 1L
          its_end   <- lsu$start - 1L
        }
        res[[op_i]] <- data.table(
          scaffold  = x$scaffold[1],
          operon_id = paste0("operon_", op_i),
          strand    = st,
          operon_start = min(pending_ssu$start, lsu$start),
          operon_end   = max(pending_ssu$end,   lsu$end),
          ssu_start = pending_ssu$start, ssu_end = pending_ssu$end,
          its_start = its_start,         its_end = its_end,
          lsu_start = lsu$start,         lsu_end = lsu$end)
        pending_ssu <- NULL
      }
    }
  }
  rbindlist(res)
}

## Extract a region (reverse-complement if on the minus strand)
get_region <- function(seqs, ops, region, start_col, end_col) {
  ok <- !is.na(ops[[start_col]]) & ops[[start_col]] <= ops[[end_col]]
  ops <- ops[ok]
  if (nrow(ops) == 0) { return(DNAStringSet()) }
  s <- subseq(seqs[ops$scaffold], start = ops[[start_col]], end = ops[[end_col]])
  minus <- ops$strand == "-"
  s[minus] <- reverseComplement(s[minus])
  names(s) <- paste0(ops$scaffold, "___", ops$operon_id, "__", region)
  s
}


## Load data
seqs <- readDNAStringSet(file.path(inp_dir, "JGI_scaffolds.fasta"))
names(seqs) <- sub(" .*", "", names(seqs))

hits_ssu <- read_tblout(file.path(inp_dir, "cmsearch_SSU.tbl"))
hits_lsu <- read_tblout(file.path(inp_dir, "cmsearch_LSU.tbl"))
if (nrow(hits_ssu) == 0 || nrow(hits_lsu) == 0) { stop("No SSU or no LSU hits found") }

hits <- rbind(
  data.table(hits_ssu, feature = "SSU"),
  data.table(hits_lsu, feature = "LSU"))

## Keep non-truncated hits on scaffolds that have both SSU and LSU
hits <- hits[ trunc %in% "no" ]
both <- intersect(hits[feature == "SSU"]$scaffold, hits[feature == "LSU"]$scaffold)
hits <- hits[ scaffold %in% both ]
if (nrow(hits) == 0) { stop("No scaffolds with both SSU and LSU hits") }

hits[ , start := pmin(seq_from, seq_to) ]
hits[ , end   := pmax(seq_from, seq_to) ]

## Find operons
OPS <- rbindlist(lapply(split(hits, by = "scaffold"), find_operons))
message("Scaffolds with complete operons: ", uniqueN(OPS$scaffold), "; operons: ", nrow(OPS))

## Export
fwrite(OPS, file.path(out_dir, "cmsearch__operons.tsv"), sep = "\t")

REG <- list(
  FullOperon = get_region(seqs, OPS, "operon", "operon_start", "operon_end"),
  SSU        = get_region(seqs, OPS, "SSU",    "ssu_start",    "ssu_end"),
  ITS        = get_region(seqs, OPS, "ITS",    "its_start",    "its_end"),
  LSU        = get_region(seqs, OPS, "LSU",    "lsu_start",    "lsu_end"))

for (r in names(REG)) {
  writeXStringSet(REG[[r]], filepath = file.path(out_dir, paste0("cmsearch__", r, ".fasta")), width = 200000)
}
