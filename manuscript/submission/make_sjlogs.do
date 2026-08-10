*! version 1.1.0  05jul2026
*! Generates the Stata Journal log excerpts (.log.tex) included in the
*! manuscript's stlog environments, using the sjlog command from the sjlatex
*! package.
*!
*! Layout-aware: run from the repository root (where code/ and examples/
*! exist) or from the flat submission packet (where the ado and .dta files sit
*! in the current directory).  Either way the commands shown in the article
*! reproduce exactly as typed.
clear all
set more off
set linesize 79

capture confirm file "code/nnsreg.ado"
if _rc == 0 {
    // Repository layout: stage a flat mirror so logged paths are packet-flat.
    cap mkdir "work"
    cap mkdir "work/packet_stage"
    foreach f in nnsreg.ado nnsreg_p.ado nnsreg_estat.ado nnsreg.sthlp {
        copy "code/`f'" "work/packet_stage/`f'", replace
    }
    foreach f in texas_education texas_water clustered_peak_dip {
        copy "examples/`f'.dta" "work/packet_stage/`f'.dta", replace
    }
    local outdir "manuscript"
    cd "work/packet_stage"
    adopath ++ "."
    local return_to "../.."
}
else {
    // Flat submission-packet layout: run in place.
    adopath ++ "."
    local outdir "."
    local return_to "."
}

which nnsreg

// --- Excerpt 1: education example, default connect fit and model comparison
sjlog using nnsreg_ex_edu, replace
use texas_education, clear
nnsreg staar_pass spend_pupil, order(3) method(connect) noplot
estimates store edu_connect
quietly nnsreg staar_pass spend_pupil, order(3) method(ols) noplot
estimates store edu_ols
estimates table edu_connect edu_ols, b(%9.3f) se(%9.3f) stats(r2 rmse N)
sjlog close, replace

// --- Excerpt 2: education plateau test via lincom on the stored spline
sjlog using nnsreg_ex_plateau, replace
estimates restore edu_ols
lincom spend_pupil + knot1 + knot2 + knot3 + knot4 + knot5 + knot6 ///
    + knot7 + knot8 + knot9
sjlog close, replace

// --- Excerpt 3: water example, slope contrast across the threshold
sjlog using nnsreg_ex_water, replace
use texas_water, clear
quietly nnsreg water_loss pipe_age, order(3) method(ols) noplot
lincom pipe_age
lincom pipe_age + knot1 + knot2 + knot3 + knot4 + knot5 + knot6 ///
    + knot7 + knot8 + knot9
sjlog close, replace

// --- Excerpt 4: the slopes option reports every segment at once
sjlog using nnsreg_ex_slopes, replace
use texas_education, clear
nnsreg staar_pass spend_pupil, order(2) method(ols) noplot slopes
sjlog close, replace

if "`return_to'" != "." {
    cd "`return_to'"
    foreach f in nnsreg_ex_edu nnsreg_ex_plateau nnsreg_ex_water nnsreg_ex_slopes {
        copy "work/packet_stage/`f'.log.tex" "`outdir'/`f'.log.tex", replace
    }
}

di "sjlog excerpts written."
