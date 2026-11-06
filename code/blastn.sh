#!/usr/bin/env bash
#
# BLASTn search of query sequences (SSU, ITS or LSU regions) against the EUKARYOME reference database
#
# Usage:
#   ./run_blastn.sh <query.fasta> <EUKARYOME.fasta> <output_prefix> [threads]
#
# Arguments:
#   query.fasta      - query sequences (e.g., ITS, SSU or LSU part of representative OTU sequences)
#   EUKARYOME.fasta  - reference sequences; headers are expected in the semicolon-delimited EUKARYOME format
#   output_prefix    - prefix for the output file (<output_prefix>.m8.gz)
#   threads          - number of CPU threads (default: 4)
#
# Output:
#   Tab-delimited BLAST table (outfmt 6) with up to 10 hits per query, columns:
#   qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
#
# Requirements: BLAST+ >= v2.16.0

set -euo pipefail

QUERY="$1"
REF="$2"
OUT="$3"
THREADS="${4:-4}"

## Build the BLAST database from the reference FASTA (skipped if it already exists)
DB="${REF%.*}_blastdb"
if [[ ! -f "${DB}.nsq" && ! -f "${DB}.00.nsq" ]]; then
  makeblastdb \
    -in "$REF" \
    -dbtype nucl \
    -out "$DB"
fi

## BLASTn search
##   -task blastn         - traditional BLASTn (more sensitive for divergent sequences than the default megablast)
##   -strand both         - search both strands
##   -max_target_seqs 10  - keep up to 10 reference hits per query
##   -max_hsps 1          - keep a single (best) HSP per query-reference pair
##   All other parameters are BLAST+ defaults
##   (word_size 11, reward 2, penalty -3, gapopen 5, gapextend 2, evalue 10, DUST filtering on)
blastn \
  -task blastn \
  -query "$QUERY" \
  -db "$DB" \
  -strand both \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$THREADS" \
  -out "${OUT}.m8"

## Compress results
gzip -f "${OUT}.m8"

echo "Done: ${OUT}.m8.gz"
