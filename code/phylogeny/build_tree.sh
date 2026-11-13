#!/usr/bin/env bash
#
# Alignment, trimming and maximum-likelihood tree inference for SSU-ITS-LSU sequences
#
# Usage:
#   ./build_tree.sh <sequences.fasta> <output_prefix> [threads]
#
# Arguments:
#   sequences.fasta - unaligned sequences (novel lineage reads + reference taxa + outgroup)
#   output_prefix   - prefix for output files
#   threads         - number of CPU threads (default: 4)
#
# Steps:
#   1. MAFFT alignment
#   2. Manual curation in AliView v1.26 (not scripted): trimming of overhanging ends,
#      removal of introns, ITS regions and the 5.8S rRNA gene, correction of misalignments.
#      Save the curated alignment as <output_prefix>.curated.fasta and re-run the script;
#      if the file is absent, the raw MAFFT alignment is used
#   3. ClipKIT trimming of phylogenetically uninformative sites
#   4. IQ-TREE with ModelFinder, 1000 ultrafast bootstrap and 1000 SH-aLRT replicates
#
# Output:
#   <output_prefix>.mafft.fasta               - raw alignment
#   <output_prefix>.clipkit.fasta             - trimmed alignment
#   <output_prefix>.clipkit.fasta.treefile    - ML tree (Newick) with SH-aLRT/UFBoot support
#   <output_prefix>.clipkit.fasta.iqtree      - IQ-TREE report (selected model, log-likelihood)
#
# Requirements: MAFFT v7.526, ClipKIT v1.4.0, IQ-TREE v2.2.5

set -euo pipefail

INPUT="$1"
PREFIX="$2"
THREADS="${3:-4}"

## 1. Alignment
if [[ ! -f "${PREFIX}.mafft.fasta" ]]; then
  mafft --auto --thread "$THREADS" "$INPUT" > "${PREFIX}.mafft.fasta"
fi

## 2. Manually curated alignment (if available)
ALN="${PREFIX}.mafft.fasta"
if [[ -f "${PREFIX}.curated.fasta" ]]; then
  ALN="${PREFIX}.curated.fasta"
fi
echo "Alignment used: ${ALN}"

## 3. Trim sites with a gap proportion above 0.85
clipkit "$ALN" -m gappy -g 0.85 -o "${PREFIX}.clipkit.fasta"

## 4. Tree inference
##   -m TEST      - model selection with ModelFinder
##   -bb 1000     - ultrafast bootstrap replicates
##   -alrt 1000   - SH-like approximate likelihood ratio test replicates
iqtree \
  -s "${PREFIX}.clipkit.fasta" \
  -m TEST \
  -bb 1000 -alrt 1000 \
  -T "$THREADS"

echo "Done: ${PREFIX}.clipkit.fasta.treefile"
