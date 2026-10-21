#!/usr/bin/env bash
# Build every per-locus genomic-context figure with ALL contexts kept, in the style settled
# on for st492_synteny_three_allctx_2col, and verify each one against the table it was drawn
# from before declaring success.
#
#   1. build    lift the five-variant cap, apply the audited corrections in memory and write
#               st492_synteny_variants_allctx8.csv (one row per arrow),
#               st492_synteny_contexts_allctx8.csv (one row per gene copy) and
#               st492_synteny_dropped_allctx8.csv (hits counted once as duplicate contigs)
#   2. plot     one figure per locus (and, with ST492_NOCOUNTRIES=1, the set without the
#               country breakdowns), plus st492_synteny_<locus>_allctx_layout.csv mapping each
#               display number to its CSV variant and st492_synteny_<locus>_allctx_caption.txt,
#               the caption text, written fresh for every locus rendered
#   3. verify   check every rendered file against the tables, its own scale bar and page size;
#               --upstream also recounts copies from marker_hits and the Bakta assemblies
#
# The tables go to ST492_TABLEDIR (default: the figure folder) and every step reads them from
# there, so a scratch run (ST492_OUTDIR and ST492_TABLEDIR pointing elsewhere) never rewrites
# the shared tables. Each rendered locus gets a FRESH st492_synteny_<locus>_allctx_caption.txt:
# the old one is deleted before the render, so a failed render cannot leave stale text behind.
#
# Nothing here touches the capped builder, the per-locus capped figures, or the combined
# figures. The builder is read and patched in memory; every output carries an _allctx name.
#
# Usage:  bash src/scripts/run_st492_synteny_allctx.sh [--skip-build] [--dpi N] [--loci A,B]
#   --loci renders and verifies only those loci (default: all eight). The build always
#   rewrites the tables for all eight, so figures not re-rendered are stale against them.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${ST492_RSCRIPT:-/home/bals/miniforge3/envs/r-ggtree/bin/Rscript}"
FIGDIR="${ST492_OUTDIR:-/Data/ST492/figures2026}"
TABLEDIR="${ST492_TABLEDIR:-$FIGDIR}"
VCSV="$TABLEDIR/st492_synteny_variants_allctx8.csv"
CCSV="$TABLEDIR/st492_synteny_contexts_allctx8.csv"
DPI="${ST492_DPI:-1200}"
SKIP_BUILD=0
LOCI_ARG="${ST492_LOCI:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-build) SKIP_BUILD=1; shift ;;
    --dpi) DPI="$2"; shift 2 ;;
    --loci) LOCI_ARG="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

# System R 4.6 cannot load the compiled packages these scripts need, so the conda R is
# required. Fail loudly rather than half way through a render.
[[ -x "$RSCRIPT" ]] || { echo "no Rscript at $RSCRIPT (set ST492_RSCRIPT)" >&2; exit 1; }
command -v pdftocairo >/dev/null || { echo "pdftocairo is needed to verify the drawn scale" >&2; exit 1; }

if [[ "$SKIP_BUILD" -eq 0 ]]; then
  echo "== 1/3  building the all-contexts tables =="
  python3 "$HERE/build_st492_synteny_contexts_allctx.py" --outdir "$TABLEDIR"
else
  echo "== 1/3  skipped, reusing $VCSV =="
fi

echo
echo "== 2/3  rendering at ${DPI} dpi =="
ALL_LOCI="blaOXA-72,blaNDM-1,intI1,blaOXA-66,blaOXA-148,blaOXA-23,carO,oprD"
IFS=, read -r -a RLOCI <<< "${LOCI_ARG:-$ALL_LOCI}"
for g in "${RLOCI[@]}"; do rm -f "$FIGDIR/st492_synteny_${g}_allctx_caption.txt"; done
ST492_LOCI="$LOCI_ARG" ST492_DPI="$DPI" ST492_OUTDIR="$FIGDIR" ST492_ALLCTX_CSV="$VCSV" \
  ST492_ALLCTX_COPIES="$CCSV" ST492_CAPTION_TXT="" "$RSCRIPT" "$HERE/plot_st492_synteny_allctx_perlocus.R"

echo
echo "== 3/3  verifying both sets against the source table =="
rc=0
VLOCI=()
[[ -n "$LOCI_ARG" ]] && VLOCI=(--loci "$LOCI_ARG")
python3 "$HERE/verify_st492_synteny_allctx.py" --dir "$FIGDIR" --csv "$VCSV" --copies "$CCSV" --dpi "$DPI" \
  --upstream "${VLOCI[@]}" || rc=1
if [[ "${ST492_NOCOUNTRIES:-0}" == "1" ]]; then
  echo
  python3 "$HERE/verify_st492_synteny_allctx.py" --dir "$FIGDIR" --csv "$VCSV" --copies "$CCSV" --dpi "$DPI" \
    --suffix _nocountries "${VLOCI[@]}" || rc=1
fi

if [[ "$rc" -ne 0 ]]; then
  echo
  echo "VERIFICATION FAILED: the figures above do not match the table they were drawn from." >&2
  exit 1
fi
echo
echo "done: ${LOCI_ARG:-all 8 loci} in $FIGDIR (set ST492_NOCOUNTRIES=1 to also build the set without country breakdowns)"
