*! version 1.0.0  17jul2026
*! Real-data illustration for the nnsreg Stata Journal manuscript.
*!
*! Data: the motorcycle crash-helmet impact experiment (accelerometer readings
*! against time after impact).  Originally Silverman (1985, JRSS-B 47: 1-52),
*! reprinted in Fan and Gijbels (1996), and shipped by Stata Press as
*! motorcycle.dta (webuse motorcycle); it is the dataset in the [R] lpoly
*! examples.  The regressor (time) is heavily tied and clustered (94 distinct
*! values in 133 observations), the mean response has a flat pre-impact run
*! followed by a sharp trough and rebound, and the noise is strongly
*! heteroskedastic, so this is nnsreg's favorable case and a natural showcase
*! for noise(median).
*!
*! Produces:
*!   - manuscript/fig_moto_fits.pdf/.png : four fits, one method per panel
*!   - manuscript/fig_moto_knots.pdf/.png: NNS OLS fit with e(knots) overlaid
*!   - manuscript/table_moto_diag.csv    : in-sample diagnostics
*!   - manuscript/table_moto_cv.csv      : repeated k-fold held-out RMSE
*!
*! Because the mean curve is unknown, method accuracy is judged by repeated
*! k-fold cross-validation (held-out RMSE), the real-data analogue of the
*! true-curve RMSE used with the simulated designs.  npregress kernel and
*! npregress series do not score observations outside their estimation sample
*! with predict, so their held-out predictions are obtained with margins at
*! the held-out regressor values; nnsreg, regress, and mkspline+regress score
*! held-out rows directly with predict.
clear all
set more off
set scheme sj
adopath + "code"
cap mkdir "manuscript"

// Featured partition depth for the single-run figures and slope table.
local ORD = 3
// Cross-validation settings.
local CVREPS = 20
local KFOLD  = 5
local CVSEED = 90210

which nnsreg

// Grayscale-safe brand grays (print and screen): dark navy fit, mid-gray dots.
local CDOT  "gs9"
local CFIT  "27 45 85"
local CFIT2 "108 122 141"

// ---------------------------------------------------------------------------
// Load and describe
// ---------------------------------------------------------------------------
webuse motorcycle, clear
count
local Nmoto = r(N)
quietly levelsof time, local(utimes)
di as text "N = `Nmoto', distinct time values = " `: word count `utimes''

// ---------------------------------------------------------------------------
// Helper: in-sample RMSE against the observed outcome (no true curve here)
// ---------------------------------------------------------------------------
capture program drop moto_rmse
program define moto_rmse, rclass
    syntax, Fitvar(name) Yvar(name)
    tempvar e2
    quietly gen double `e2' = (`yvar' - `fitvar')^2
    quietly summarize `e2', meanonly
    return scalar rmse = sqrt(r(mean))
    quietly correlate `yvar' `fitvar'
    return scalar r2 = r(rho)^2
end

// Panel: gray dots + one dark fitted curve, one method per panel.
capture program drop moto_panel
program define moto_panel
    syntax, Fitvar(name) PTitle(string asis) PSub(string asis) GName(name)
    twoway ///
        (scatter accel time, mcolor("108 122 141%45") msymbol(Oh) msize(small)) ///
        (line `fitvar' time, sort lcolor("27 45 85") lwidth(medthick)), ///
        title(`ptitle', size(medsmall) position(11) span) ///
        subtitle(`psub', size(small) position(11) span color("108 122 141")) ///
        ytitle("") xtitle("") ///
        ylabel(-150(50)100, labsize(small) angle(horizontal)) ///
        xlabel(0(20)60, labsize(small)) ///
        legend(off) graphregion(color(white)) name(`gname', replace) nodraw
end

tempname diag
postfile `diag' str30 model double(r2 rmse_obs) nknots using ///
    "examples/moto_diag.dta", replace

// ---------------------------------------------------------------------------
// Single-run fits (in sample) and the fits figure
// ---------------------------------------------------------------------------
regress accel time
predict double fit_lin, xb
moto_rmse, fitvar(fit_lin) yvar(accel)
post `diag' ("Linear regression") (r(r2)) (r(rmse)) (0)
local rl : di %5.2f r(rmse)

quietly npregress kernel accel time
predict double fit_np, mean
moto_rmse, fitvar(fit_np) yvar(accel)
post `diag' ("npregress kernel") (r(r2)) (r(rmse)) (.)
local rk : di %5.2f r(rmse)

quietly npregress series accel time
predict double fit_series, mean
moto_rmse, fitvar(fit_series) yvar(accel)
post `diag' ("npregress series") (r(r2)) (r(rmse)) (.)

nnsreg accel time, order(`ORD') method(connect) noplot
local k_conn : word count `e(knots)'
predict double fit_conn
moto_rmse, fitvar(fit_conn) yvar(accel)
post `diag' ("NNS connect, order(`ORD')") (r(r2)) (r(rmse)) (`k_conn')
local rc : di %5.2f r(rmse)

nnsreg accel time, order(`ORD') method(ols) noplot
local k_ols : word count `e(knots)'
predict double fit_ols
moto_rmse, fitvar(fit_ols) yvar(accel)
post `diag' ("NNS OLS, order(`ORD')") (r(r2)) (r(rmse)) (`k_ols')
local ro : di %5.2f r(rmse)
estimates store moto_ols

nnsreg accel time, order(`ORD') method(ols) noise(median) noplot
predict double fit_ols_med
moto_rmse, fitvar(fit_ols_med) yvar(accel)
post `diag' ("NNS OLS, order(`ORD') median") (r(r2)) (r(rmse)) (`: word count `e(knots)'')
postclose `diag'

moto_panel, fitvar(fit_lin) ptitle("A. Linear regression") ///
    psub("In-sample RMSE = `rl'") gname(m_p1)
moto_panel, fitvar(fit_np) ptitle("B. npregress kernel") ///
    psub("In-sample RMSE = `rk'") gname(m_p2)
moto_panel, fitvar(fit_conn) ptitle("C. NNS connect, order(`ORD')") ///
    psub("In-sample RMSE = `rc'") gname(m_p3)
moto_panel, fitvar(fit_ols) ptitle("D. NNS OLS, order(`ORD')") ///
    psub("In-sample RMSE = `ro'") gname(m_p4)
graph combine m_p1 m_p2 m_p3 m_p4, cols(2) imargin(small) ///
    l1title("Acceleration (g)", size(small)) ///
    b1title("Time after impact (milliseconds)", size(small)) ///
    note("Gray dots: 133 accelerometer readings.  Dark line: fitted values.", size(vsmall)) ///
    graphregion(color(white)) xsize(7.5) ysize(6) name(moto_fits, replace)
graph export "manuscript/fig_moto_fits.pdf", replace
graph export "manuscript/fig_moto_fits.png", replace width(2200)

// Knots figure: the adaptive knots concentrate where the curve bends fastest.
estimates restore moto_ols
local klist "`e(knots)'"
twoway ///
    (scatter accel time, mcolor("108 122 141%45") msymbol(Oh) msize(small)) ///
    (line fit_ols time, sort lcolor("27 45 85") lwidth(medthick)), ///
    xline(`klist', lcolor("108 122 141") lpattern(shortdash) lwidth(vthin)) ///
    ytitle("Acceleration (g)", size(small)) ///
    xtitle("Time after impact (milliseconds)", size(small)) ///
    ylabel(-150(50)100, labsize(small) angle(horizontal)) xlabel(0(20)60, labsize(small)) ///
    legend(order(2 "NNS OLS fit, order(`ORD')") rows(1) size(small) position(5) ///
        region(lcolor(none))) ///
    note("Vertical dashed lines: the `k_ols' interior knots from e(knots); they concentrate in the trough where the curve bends fastest.", ///
        size(vsmall)) ///
    graphregion(color(white)) xsize(7) ysize(4.5) name(moto_knots, replace)
graph export "manuscript/fig_moto_knots.pdf", replace
graph export "manuscript/fig_moto_knots.png", replace width(2000)

// Diagnostics table
use "examples/moto_diag.dta", clear
format r2 rmse_obs %6.3f
list, noobs
export delimited using "manuscript/table_moto_diag.csv", replace

// Featured slope table (for a manuscript log excerpt / reference)
di _n "=== per-segment slopes on the featured NNS OLS fit ==="
webuse motorcycle, clear
nnsreg accel time, order(`ORD') method(ols) noplot slopes

// ---------------------------------------------------------------------------
// Repeated k-fold cross-validation: held-out RMSE by method
// ---------------------------------------------------------------------------
di _n "=== Repeated `KFOLD'-fold CV, `CVREPS' repeats ==="
webuse motorcycle, clear
// Fixed restricted-cubic-spline basis at sample percentile knots (the
// standard mkspline usage; knots placed by rank, not by where the curve bends).
mkspline msp = time, cubic nknots(5)

tempname cv
postfile `cv' str24 method double rmse int rep using "examples/moto_cv.dta", replace

// npregress series is excluded from the held-out comparison: its predict is
// in-sample only, and margins-based out-of-sample prediction refuses points
// outside the fitted B-spline support (returns "prediction outside valid
// range"), so it cannot be scored on the same held-out set as the other arms.
// Its in-sample fit is in table_moto_diag.csv.  npregress kernel scores
// held-out points cleanly via margins.
local methods "linear mksp kernel conn ols olsmed"
local npfail = 0

// The observations at the smallest and largest time are always kept in the
// training set (fold 0), so every held-out point is interior to the training
// range.  No smoother can honestly extrapolate beyond the data, and this keeps
// npregress (whose margins refuses out-of-range prediction) on the same
// held-out set as the arms that predict out of sample directly.
quietly summarize time, meanonly
local tmin = r(min)
local tmax = r(max)

forvalues rep = 1/`CVREPS' {
    set seed `=`CVSEED' + `rep''
    quietly gen double u = runiform()
    quietly xtile fold = u, nq(`KFOLD')
    quietly replace fold = 0 if time == `tmin' | time == `tmax'
    // One OOS prediction column per method
    foreach m of local methods {
        quietly gen double p_`m' = .
    }
    forvalues k = 1/`KFOLD' {
        // ----- arms that predict out-of-sample directly -----
        quietly regress accel time if fold != `k'
        quietly predict double _tmp, xb
        quietly replace p_linear = _tmp if fold == `k'
        drop _tmp

        quietly regress accel msp* if fold != `k'
        quietly predict double _tmp, xb
        quietly replace p_mksp = _tmp if fold == `k'
        drop _tmp

        capture quietly nnsreg accel time if fold != `k', order(`ORD') method(connect) noplot
        if _rc == 0 {
            quietly predict double _tmp
            quietly replace p_conn = _tmp if fold == `k'
            drop _tmp
        }
        capture quietly nnsreg accel time if fold != `k', order(`ORD') method(ols) noplot
        if _rc == 0 {
            quietly predict double _tmp
            quietly replace p_ols = _tmp if fold == `k'
            drop _tmp
        }
        capture quietly nnsreg accel time if fold != `k', order(`ORD') method(ols) noise(median) noplot
        if _rc == 0 {
            quietly predict double _tmp
            quietly replace p_olsmed = _tmp if fold == `k'
            drop _tmp
        }

        // ----- npregress arms: held-out predictions via margins at() -----
        levelsof time if fold == `k', local(hv)
        capture {
            quietly npregress kernel accel time if fold != `k'
            quietly margins, at(time=(`hv')) predict(mean)
            matrix P = r(b)
            local i = 0
            foreach v of local hv {
                local ++i
                quietly replace p_kernel = P[1, `i'] if fold == `k' & abs(time - `v') < 0.001
            }
        }
        if _rc local ++npfail
    }
    // One held-out RMSE per method for this repeat (over all N rows)
    foreach m of local methods {
        tempvar e2
        quietly gen double `e2' = (accel - p_`m')^2
        quietly summarize `e2', meanonly
        post `cv' ("`m'") (sqrt(r(mean))) (`rep')
        drop `e2' p_`m'
    }
    drop u fold
    if mod(`rep', 5) == 0 di as text "  CV repeat `rep' of `CVREPS'"
}
postclose `cv'
di as text "npregress margins failures across all CV folds: `npfail'"
assert `npfail' == 0

use "examples/moto_cv.dta", clear
gen msort = .
local i = 0
foreach m of local methods {
    local ++i
    quietly replace msort = `i' if method == "`m'"
}
gen str30 label = ""
replace label = "Linear regression"          if method == "linear"
replace label = "Restricted cubic spline (5)" if method == "mksp"
replace label = "npregress kernel"            if method == "kernel"
replace label = "NNS connect, order(`ORD')"     if method == "conn"
replace label = "NNS OLS, order(`ORD')"         if method == "ols"
replace label = "NNS OLS, order(`ORD') median"  if method == "olsmed"

collapse (mean) mean_cv_rmse = rmse (sd) sd_cv_rmse = rmse (count) reps = rmse, ///
    by(method label msort)
sort msort
format mean_cv_rmse sd_cv_rmse %6.3f
list label mean_cv_rmse sd_cv_rmse reps, noobs
export delimited label mean_cv_rmse sd_cv_rmse reps using ///
    "manuscript/table_moto_cv.csv", replace

di _n "MOTORCYCLE EXAMPLE COMPLETE"
