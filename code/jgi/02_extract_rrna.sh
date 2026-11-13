#!/usr/bin/env bash
#
# Combine downloaded IMG/MER scaffolds and locate SSU and LSU rRNA genes with Infernal
#
# Usage:
#   ./02_extract_rrna.sh <downloads_dir> <output_dir> [threads]
#
# Arguments:
#   downloads_dir - directory with *.tar.gz archives from the IMG/MER "Export Scaffold Data" job
#   output_dir    - output directory
#   threads       - number of CPU threads for cmsearch (default: 4)
#
# Output:
#   <output_dir>/JGI_scaffolds.fasta - all scaffolds; spaces in headers are replaced with "___"
#                                      (e.g., >3300020163___assembled___Ga0194039_1011500)
#   <output_dir>/RF01960.cm, RF02543.cm - Rfam covariance models (eukaryotic SSU and LSU rRNA)
#   <output_dir>/cmsearch_SSU.tbl, cmsearch_LSU.tbl - cmsearch hit tables (both strands)
#
# Next step: Rscript 03_rrna_operons.R <output_dir> <output_dir>
#
# Requirements: seqkit, GNU parallel, Infernal v1.1.5, wget

set -euo pipefail

DLDIR="$1"
OUTDIR="$2"
THREADS="${3:-4}"

mkdir -p "$OUTDIR"

## Unpack scaffold sequences from all archives
find "$DLDIR" -type f -name "*.tar.gz" | sort \
  | parallel -j 1 "tar -xzOf {} --wildcards '*/fasta.fna'" \
  | seqkit seq -w 0 \
  | sed '/^>/ s/ /___/g' \
  > "$OUTDIR"/JGI_scaffolds.fasta

seqkit stats "$OUTDIR"/JGI_scaffolds.fasta

## Rfam covariance models
##   RF01960 - eukaryotic SSU rRNA
##   RF02543 - eukaryotic LSU rRNA
for CM in RF01960 RF02543; do
  if [[ ! -f "$OUTDIR"/${CM}.cm ]]; then
    wget -q -O "$OUTDIR"/${CM}.cm "https://rfam.org/family/${CM}/cm"
  fi
done

## Search both strands; keep hits in tabular format
cmsearch --cpu "$THREADS" --noali -o /dev/null \
  --tblout "$OUTDIR"/cmsearch_SSU.tbl \
  "$OUTDIR"/RF01960.cm "$OUTDIR"/JGI_scaffolds.fasta

cmsearch --cpu "$THREADS" --noali -o /dev/null \
  --tblout "$OUTDIR"/cmsearch_LSU.tbl \
  "$OUTDIR"/RF02543.cm "$OUTDIR"/JGI_scaffolds.fasta

echo "Done: ${OUTDIR}"
