#!/usr/bin/env bash
#
# BLASTn of scaffold rRNA regions and of their 5' and 3' ends against EUKARYOME
# (taxonomic assignment and detection of chimeric scaffolds)
#
# Usage:
#   ./04_blast_ends.sh <regions_dir> <EUKARYOME_SSU.fasta> <EUKARYOME_LSU.fasta> <EUKARYOME_ITS.fasta> <output_dir> [threads]
#
# Arguments:
#   regions_dir  - directory with cmsearch__{SSU,ITS,LSU,FullOperon}.fasta (from 03_rrna_operons.R)
#   EUKARYOME_*  - EUKARYOME reference sequences for SSU, LSU and ITS (https://eukaryome.org)
#   output_dir   - output directory
#   threads      - number of CPU threads (default: 4)
#
# Queries:
#   SSU, LSU             - complete regions
#   SSU/LSU first 500 bp - 5' end of the region
#   SSU/LSU last 500 bp  - 3' end of the region
#   ITS, full operon     - against the ITS reference
#
# Output: <output_dir>/<query>.m8.gz (BLAST tabular format; see ../blastn.sh)
#
# Interpretation (as in the manuscript):
#   - scaffolds whose parts match taxonomically incongruent groups at kingdom or phylum level are considered chimeric
#   - phylum-level identification requires >= 94.0% identity for SSU and >= 92.0% for LSU
#   - species-level identity requires a 100% identical match over at least 500 bases
#
# Requirements: seqkit, BLAST+ v2.16.0

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLASTN="${SCRIPT_DIR}/../blastn.sh"

IN="$1"
REF_SSU="$2"
REF_LSU="$3"
REF_ITS="$4"
OUT="$5"
THREADS="${6:-4}"

mkdir -p "$OUT"/queries

## First and last 500 bases of SSU and LSU
for R in SSU LSU; do
  cp "$IN"/cmsearch__${R}.fasta "$OUT"/queries/${R}.fasta
  seqkit subseq -r 1:500    -w 0 "$IN"/cmsearch__${R}.fasta > "$OUT"/queries/${R}_first_500.fasta
  seqkit subseq -r -500:-1  -w 0 "$IN"/cmsearch__${R}.fasta > "$OUT"/queries/${R}_last_500.fasta
done
cp "$IN"/cmsearch__ITS.fasta        "$OUT"/queries/ITS.fasta
cp "$IN"/cmsearch__FullOperon.fasta "$OUT"/queries/FullOperon.fasta

## BLASTn
for Q in SSU SSU_first_500 SSU_last_500; do
  "$BLASTN" "$OUT"/queries/${Q}.fasta "$REF_SSU" "$OUT"/${Q} "$THREADS"
done

for Q in LSU LSU_first_500 LSU_last_500; do
  "$BLASTN" "$OUT"/queries/${Q}.fasta "$REF_LSU" "$OUT"/${Q} "$THREADS"
done

for Q in ITS FullOperon; do
  "$BLASTN" "$OUT"/queries/${Q}.fasta "$REF_ITS" "$OUT"/${Q} "$THREADS"
done

echo "Done: ${OUT}"
