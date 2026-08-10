#!/bin/bash
# Rebuilds the flat Stata Journal submission packet in manuscript/submission/
# from the canonical repository files (the manuscript/ folder is the single
# source of truth).
#
# The packet is lean and self-contained: nnsreg.tex compiles in place, and the
# path-flattened do-files reproduce every example from inside the directory.
# It holds only what a Stata Journal submission needs -- manuscript sources and
# style files, the compiled PDF, the package ado/help files, the reproduction
# do-files and datasets, the output excerpts and figures the manuscript uses,
# a readme, and the cover letter. It contains no LaTeX build artifacts, no PNG
# figure twins, no intermediate CSVs, and no regenerable intermediate datasets.
#
# Usage: ./sync_submission.sh
set -euo pipefail
cd "$(dirname "$0")"

DEST="manuscript/submission"

flatten() {
    # $1 = source do-file, $2 = destination; rewrite repo paths to flat paths.
    sed -e 's/adopath + "code"/adopath ++ "."/' \
        -e '/^cap mkdir "manuscript"$/d' \
        -e '/^cap mkdir "examples"$/d' \
        -e 's|"examples/|"|g' \
        -e 's|"manuscript/|"|g' \
        "$1" > "$2"
}

# Start from an empty packet so no stale or orphaned file can survive.
rm -rf "$DEST"
mkdir -p "$DEST"

# --- Manuscript source (single self-contained file) and SJ LaTeX files ---
cp manuscript/nnsreg.tex manuscript/nnsreg.bib "$DEST/"
cp manuscript/statapress.cls manuscript/sj.sty manuscript/stata.sty \
   manuscript/pagedims.sty manuscript/sj.bst "$DEST/"

# --- Stata output excerpts included by the manuscript ---
cp manuscript/nnsreg_ex_edu.log.tex manuscript/nnsreg_ex_plateau.log.tex \
   manuscript/nnsreg_ex_water.log.tex manuscript/nnsreg_ex_slopes.log.tex "$DEST/"

# --- Figures referenced by the manuscript (PDF only; the Stata Journal prefers
#     PDF and grayscale). Stata's macOS PDF export references Helvetica without
#     embedding it, so embed a substitute font via ghostscript (idempotent). ---
for f in fig_edu_fits fig_edu_coef fig_edu_sensitivity \
         fig_water_fits fig_water_knots fig_water_sensitivity \
         fig_cluster_fits fig_cluster_sensitivity fig_mc_rmse \
         fig_moto_fits fig_moto_knots; do
    if command -v gs > /dev/null && ! pdffonts "manuscript/$f.pdf" | grep -q "yes yes"; then
        gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -dEmbedAllFonts=true \
           -dSubsetFonts=true -dCompatibilityLevel=1.5 \
           -o "manuscript/$f.embedded.pdf" \
           -c "<</NeverEmbed []>> setdistillerparams" -f "manuscript/$f.pdf" \
           && mv "manuscript/$f.embedded.pdf" "manuscript/$f.pdf"
    fi
    cp "manuscript/$f.pdf" "$DEST/"
done

# --- Compiled manuscript ---
cp manuscript/nnsreg.pdf "$DEST/"

# --- Package files (command, predict, estat dispatcher, help) ---
cp code/nnsreg.ado code/nnsreg_p.ado code/nnsreg_estat.ado code/nnsreg.sthlp "$DEST/"

# --- Reproducibility do-files, path-flattened for the flat packet.
#     make_sjlogs.do detects flat-vs-repo layout at runtime, so it is copied
#     verbatim. ---
flatten examples/simulated_data.do "$DEST/simulated_data.do"
flatten examples/test_nnsreg.do "$DEST/test_nnsreg.do"
flatten examples/motorcycle_example.do "$DEST/motorcycle_example.do"
flatten examples/robustness.do "$DEST/robustness.do"
cp examples/make_sjlogs.do "$DEST/make_sjlogs.do"

# --- Example datasets used by the do-files (regenerable by simulated_data.do) ---
cp examples/texas_education.dta examples/texas_water.dta \
   examples/clustered_peak_dip.dta "$DEST/"

# --- Readme and cover letter (source lives in manuscript/) ---
cp manuscript/readme.txt "$DEST/readme.txt"
cp manuscript/cover_letter.txt "$DEST/cover_letter.txt"

echo "Lean submission packet rebuilt in $DEST"
echo "Files: $(ls -1 "$DEST" | wc -l | tr -d ' ')"
