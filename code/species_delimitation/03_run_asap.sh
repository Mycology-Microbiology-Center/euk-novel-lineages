#!/usr/bin/env bash
#
# Run ASAP (Assemble Species by Automatic Partitioning; Puillandre et al. 2021)
# using iTaxoTools ASAPy v0.1.2 (https://github.com/iTaxoTools/ASAPy)
#
# Usage:
#   ./03_run_asap.sh <input> [output_dir]
#
#   input      - distance matrix in MEGA CSV format (from 02_distances_to_mega.R; set ASAP_GENERAL_MEGA=true)
#                or an aligned FASTA file (gaps allowed, all sequences of the same width)
#   output_dir - default: ./asap_results/<input_name>
#
# Example (all lineages, as used in the study):
#   find mega/ -name "*.csv" | parallel -j1 \
#     "ASAP_GENERAL_SEQUENCE_LEN=600 ASAP_GENERAL_MEGA=true ./03_run_asap.sh {} asap/{/.}"
#
# ASAP parameters (environment variables; defaults match the ASAPy GUI):
#   ASAP_GENERAL_ALL           true|false  Generate all partition files (default: true)
#   ASAP_GENERAL_MEGA          true|false  Input is a MEGA distance CSV, not FASTA (default: false)
#   ASAP_GENERAL_SEQUENCE_LEN  int         Sequence length for coalescent scaling (default: 600)
#   ASAP_ADVANCED_NUMBER       int         Number of top-ranked partitions to report (default: 10)
#   ASAP_ADVANCED_SEUIL_PVALUE float       p-value threshold for recursive splits (default: 0.01)
#   ASAP_ADVANCED_SEED         int         RNG seed; -1 = not fixed (default: -1)
#   ASAP_DISTANCE_METHOD       int         0 = K80, 1 = JC69, 2 = TN93, 3 = p-distance (default: 1; ignored for MEGA input)
#   ASAP_DISTANCE_RATE         float       K80 transition/transversion ratio (default: 2.0)
#
# Notes:
#   - ASAP returns several candidate partitions ranked by the ASAP score.
#     The top-ranked partition is a statistical default; compare it with the phylogeny.
#   - For FASTA input, positions with a gap or N in either sequence are excluded from each pairwise distance.
#
# Installation:
#   mamba create -y -n ASAP python=3.11 compilers && mamba activate ASAP
#   git clone https://github.com/iTaxoTools/ASAPy && cd ASAPy
#   export CFLAGS="${CFLAGS:-} -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types"
#   python -m pip install . -f packages.html

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
fi

INPUT="$1"
if [[ ! -f "$INPUT" ]]; then
  echo "Error: not a file: $INPUT" >&2
  exit 1
fi

STEM="$(basename "$INPUT")"
STEM="${STEM%.*}"
OUTDIR="${2:-./asap_results/${STEM}}"

: "${ASAP_GENERAL_ALL:=true}"
: "${ASAP_GENERAL_MEGA:=false}"
: "${ASAP_GENERAL_SEQUENCE_LEN:=600}"
: "${ASAP_ADVANCED_NUMBER:=10}"
: "${ASAP_ADVANCED_SEUIL_PVALUE:=0.01}"
: "${ASAP_ADVANCED_SEED:=-1}"
: "${ASAP_DISTANCE_METHOD:=1}"
: "${ASAP_DISTANCE_RATE:=2.0}"

if ! python3 -c "from itaxotools.asapy.core import PartitionAnalysis" 2>/dev/null; then
  echo "Error: ASAPy is not importable, install it first." >&2
  exit 1
fi

echo "=== ASAP: ${INPUT} -> ${OUTDIR} ==="

export ASAP_RUN_INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
export ASAP_RUN_OUT="$OUTDIR"
export ASAP_GENERAL_ALL ASAP_GENERAL_MEGA ASAP_GENERAL_SEQUENCE_LEN
export ASAP_ADVANCED_NUMBER ASAP_ADVANCED_SEUIL_PVALUE ASAP_ADVANCED_SEED
export ASAP_DISTANCE_METHOD ASAP_DISTANCE_RATE

## PartitionAnalysis.launch() runs ASAP in a subprocess
## fetch() copies results from the temp dir
python3 <<'PY'
import os
import shutil
from pathlib import Path

from itaxotools.asapy.core import PartitionAnalysis

def as_bool(name):
    v = os.environ[name].strip().lower()
    if v in ("1", "true", "yes", "on"):
        return True
    if v in ("0", "false", "no", "off"):
        return False
    raise SystemExit(f"{name} must be true or false, got {os.environ[name]!r}")

outdir = Path(os.environ["ASAP_RUN_OUT"])

a = PartitionAnalysis(os.environ["ASAP_RUN_INPUT"])
a.params.general.all              = as_bool("ASAP_GENERAL_ALL")
a.params.general.mega             = as_bool("ASAP_GENERAL_MEGA")
a.params.general.sequence_length  = int(os.environ["ASAP_GENERAL_SEQUENCE_LEN"])
a.params.advanced.number          = int(os.environ["ASAP_ADVANCED_NUMBER"])
a.params.advanced.seuil_pvalue    = float(os.environ["ASAP_ADVANCED_SEUIL_PVALUE"])
a.params.advanced.seed            = int(os.environ["ASAP_ADVANCED_SEED"])
a.params.distance.method          = int(os.environ["ASAP_DISTANCE_METHOD"])
a.params.distance.rate            = float(os.environ["ASAP_DISTANCE_RATE"])

a.launch()
if outdir.exists():
    shutil.rmtree(outdir)
outdir.parent.mkdir(parents=True, exist_ok=True)
a.fetch(str(outdir))
print("Saved:", outdir)
PY
