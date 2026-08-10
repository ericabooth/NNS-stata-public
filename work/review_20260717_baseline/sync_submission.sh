#!/bin/bash
# Rebuilds the flat Stata Journal submission packet in manuscript/submission/
# from the canonical repo files. The packet do-files are path-flattened copies
# of examples/simulated_data.do and examples/test_nnsreg.do, so every example
# and figure reproduces from inside the packet directory as typed.
#
# Usage:
#   ./sync_submission.sh          sync packet from canonical files
#   ./sync_submission.sh stage D  write flattened packet into directory D only
set -euo pipefail
cd "$(dirname "$0")"

DEST="${2:-manuscript/submission}"
MODE="${1:-sync}"

flatten() {
    # $1 = source do-file, $2 = destination
    sed -e 's/adopath + "code"/adopath ++ "."/' \
        -e '/^cap mkdir "manuscript"$/d' \
        -e '/^cap mkdir "examples"$/d' \
        -e 's|"examples/|"|g' \
        -e 's|"manuscript/|"|g' \
        "$1" > "$2"
}

mkdir -p "$DEST"

# Package files
cp code/nnsreg.ado code/nnsreg_p.ado code/nnsreg.sthlp "$DEST/"

# Reproducibility scripts, path-flattened for the flat packet layout.
# simulated_data.do and test_nnsreg.do write to examples/ and manuscript/
# subpaths, so they are path-flattened.  make_sjlogs.do is layout-aware
# (it detects flat vs repo at runtime), so it is copied verbatim.
flatten examples/simulated_data.do "$DEST/simulated_data.do"
flatten examples/test_nnsreg.do "$DEST/test_nnsreg.do"
flatten examples/motorcycle_example.do "$DEST/motorcycle_example.do"
cp examples/make_sjlogs.do "$DEST/make_sjlogs.do"

# The R fidelity check (r_fidelity.R, moto_for_r.csv) is a supplementary
# cross-language alignment check kept in the development repository only, not
# shipped in the Stata Journal packet, which stays Stata-only.

# Canonical simulated datasets (regenerable via simulated_data.do)
cp examples/texas_education.dta examples/texas_water.dta \
   examples/clustered_peak_dip.dta "$DEST/"

if [ "$MODE" = "stage" ]; then
    echo "Staged flattened packet in $DEST"
    exit 0
fi

# Manuscript sources and SJ style files
cp manuscript/main.tex manuscript/nnsreg.tex manuscript/nnsreg.bib "$DEST/"
cp manuscript/statapress.cls manuscript/sj.sty manuscript/stata.sty \
   manuscript/pagedims.sty manuscript/sj.bst "$DEST/"

# Real Stata output excerpts included by the manuscript
cp manuscript/nnsreg_ex_edu.log.tex manuscript/nnsreg_ex_plateau.log.tex \
   manuscript/nnsreg_ex_water.log.tex manuscript/nnsreg_ex_slopes.log.tex "$DEST/"

# Figures referenced by the manuscript (PDF, vector) plus PNG twins.
# Stata's macOS PDF export references Helvetica without embedding it, which
# strict PDF pipelines render as missing text; run each figure through
# ghostscript to embed a substitute font (idempotent).
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
    cp "manuscript/$f.png" "$DEST/"
done

# Tables exported by the replication scripts
cp manuscript/table_model_diagnostics.csv manuscript/table_mc_summary.csv \
   manuscript/table_postest.csv "$DEST/"
cp manuscript/table_moto_diag.csv manuscript/table_moto_cv.csv "$DEST/"

# Compiled manuscript
cp manuscript/main.pdf "$DEST/"

# Remove stale artifacts from earlier drafts
rm -f "$DEST"/r_fidelity.R "$DEST"/moto_for_r.csv "$DEST"/table_moto_rfidelity.csv
rm -f "$DEST"/fig_edu_combined_compare.* "$DEST"/fig_water_combined_compare.* \
      "$DEST"/fig_cluster_combined_compare.* "$DEST"/fig_edu_coefplot.* \
      "$DEST"/fig_edu_connect_o2.* "$DEST"/fig_edu_connect_o3.* \
      "$DEST"/fig_edu_ols_o3.* "$DEST"/fig_water_connect_o2.* \
      "$DEST"/fig_water_connect_o3.* "$DEST"/fig_water_ols_o3.* \
      "$DEST"/fig_cluster_sensitivity_old.*

echo "Submission packet synced to $DEST"
