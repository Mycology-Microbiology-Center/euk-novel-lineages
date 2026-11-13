#!/usr/bin/env bash
#
# Pairwise ITS sequence similarities within each lineage (input for ASAPy and hierarchical clustering)
#
# Usage:
#   ./01_pairwise_distances.sh <alignments_dir> <output_dir> [threads]
#
# Arguments:
#   alignments_dir - directory with per-lineage ITS alignments (*.fasta or *.fasta.gz)
#   output_dir     - output directory
#   threads        - number of CPU threads for VSEARCH (default: 4)
#
# Steps (for each alignment):
#   1. Remove gaps, convert to upper case, discard sequences < 10 nt
#   2. Remove duplicated sequences; when several records share the same sequence,
#      keep the type sequence (header contains "type"), otherwise the smallest accession ID
#   3. All-vs-all global alignment with VSEARCH (every query against every other sequence)
#
# Output:
#   <output_dir>/seqs/<lineage>.fasta       - unaligned, dereplicated sequences
#   <output_dir>/dist/<lineage>.blast6.gz   - pairwise identities (BLAST tabular format)
#
# Lineage names are taken from file names, with the leading underscore and
# the threshold suffix removed (e.g., `_Anhoristes2_0.fasta.gz` -> `Anhoristes`).
#
# Requirements: seqkit, VSEARCH (>= 2.x)

set -euo pipefail

INDIR="$1"
OUTDIR="$2"
THREADS="${3:-4}"

mkdir -p "$OUTDIR"/seqs "$OUTDIR"/dist

shopt -s nullglob
for f in "$INDIR"/*.fasta "$INDIR"/*.fasta.gz; do

  name="$(basename "$f")"
  name="${name#_}"
  name="$(echo "$name" | sed -E 's/\.fasta(\.gz)?$//; s/[0-9]+_[0-9]+$//')"
  echo "== ${name}"

  ## Degap, dereplicate (type sequence first, then accession ID)
  seqkit seq --remove-gaps --upper-case --min-len 10 -w 0 "$f" \
    | seqkit fx2tab \
    | awk -F'\t' -v OFS='\t' '{
        id = $1; sub(/;$/, "", id); split(id, a, ";")
        print (id ~ /type/ ? 1 : 2), a[1], id, $2 }' \
    | LC_ALL=C sort -t$'\t' -k1,1n -k2,2 \
    | awk -F'\t' '!seen[$4]++ { print ">" $3; print $4 }' \
    > "$OUTDIR"/seqs/"$name".fasta

  ## All-vs-all global pairwise identities
  ##   --iddef 2      - identity excluding terminal gaps
  ##   --qmask dust   - mask low-complexity regions
  ##   --self         - skip self-hits
  ##   --maxaccepts 0 --maxrejects 0 - report all pairs (exhaustive search)
  vsearch \
    --usearch_global "$OUTDIR"/seqs/"$name".fasta \
    --db             "$OUTDIR"/seqs/"$name".fasta \
    --id 0 --iddef 2 --qmask dust \
    --self --maxaccepts 0 --maxrejects 0 \
    --threads "$THREADS" --quiet \
    --blast6out "$OUTDIR"/dist/"$name".blast6

  gzip -f "$OUTDIR"/dist/"$name".blast6
done

echo "Done: ${OUTDIR}/dist/"
