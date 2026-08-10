*! version 2.1.0  18jul2026
*! Replication, comparison, and testing script for the nnsreg package.
*!
*! Reproduces every figure and table in the Stata Journal manuscript:
*!   Section 1: verification tests (stored results, predict, error handling)
*!   Sections 2-4: the three simulated examples (education, water, clustered),
*!     each with a small-multiples fit-comparison figure (one method per
*!     panel, so no overlapping curves), sensitivity figures, and, for the
*!     education example, a coefficient plot and a plateau test via lincom
*!   Section 5: Monte Carlo study across `MCREPS' fresh draws of all three
*!     designs, because a comparison based on one simulated dataset is weak
*!     evidence; exports summary table and box-plot figure
*!   Section 6: single-run diagnostics table
*!
*! The metric "RMSE against the true curve" is available because the data are
*! simulated: it measures recovery of the known signal, so a method that fits
*! noise slightly better (lower observed RMSE) can still score worse on it.
clear all
macro drop _all
set more off
set scheme sj

adopath + "code"
cap mkdir "manuscript"
cap mkdir "examples"

which nnsreg
which nnsreg_p

capture which coefplot
if _rc != 0 {
    di "coefplot not found. Installing from SSC..."
    ssc install coefplot, replace
}

// Number of Monte Carlo replications (section 5). 100 runs in a few minutes.
local MCREPS = 100

// Generates canonical datasets and defines the reusable nns_simdata program.
do "examples/simulated_data.do"

// ---------------------------------------------------------------------------
// Helper programs
// ---------------------------------------------------------------------------

// calc_metrics: in-sample R2, RMSE vs observed outcome, RMSE vs true curve.
capture program drop calc_metrics
program define calc_metrics, rclass
    syntax, Fitvar(name) Yvar(name) Truevar(name)

    tempvar res_obs res_true tss
    quietly gen double `res_obs' = (`yvar' - `fitvar')^2
    quietly gen double `res_true' = (`truevar' - `fitvar')^2
    quietly summarize `res_obs', meanonly
    local rss = r(sum)
    local rmse_obs = sqrt(r(mean))
    quietly summarize `res_true', meanonly
    local rmse_true = sqrt(r(mean))
    quietly summarize `yvar', meanonly
    local ybar = r(mean)
    quietly gen double `tss' = (`yvar' - `ybar')^2
    quietly summarize `tss', meanonly

    return scalar r2 = 1 - (`rss' / r(sum))
    return scalar rmse_obs = `rmse_obs'
    return scalar rmse_true = `rmse_true'
end

// fitpanel: one small-multiples panel = gray scatter + black true curve +
// one orange fitted curve.  Keeping one method per panel avoids the
// overlapping-line clutter of a single six-layer comparison graph.
capture program drop fitpanel
program define fitpanel
    syntax, Yvar(name) Xvar(name) Truevar(name) Fitvar(name) ///
        PTitle(string asis) PSub(string asis) GName(name) [YLab(string asis)]

    twoway ///
        (scatter `yvar' `xvar', mcolor("108 122 141%35") msymbol(Oh) msize(tiny)) ///
        (line `truevar' `xvar', sort lcolor(black) lwidth(medthin)) ///
        (line `fitvar' `xvar', sort lcolor("212 69 0") lwidth(medthick)), ///
        title(`ptitle', size(medsmall) position(11) span) ///
        subtitle(`psub', size(small) position(11) span color("108 122 141")) ///
        ytitle("") xtitle("") ///
        ylabel(`ylab', labsize(small) angle(horizontal)) ///
        xlabel(, labsize(small)) ///
        legend(off) graphregion(color(white)) ///
        name(`gname', replace) nodraw
end

// Diagnostics postfile for the single-run table (manuscript table 1).
tempname diag
postfile `diag' str14 example str30 model double(r2 rmse_obs rmse_true) using ///
    "examples/model_diagnostics.dta", replace

// ---------------------------------------------------------------------------
// 1. Verification tests: stored results, predict, if/in, error handling
// ---------------------------------------------------------------------------
di _n "=== Section 1: verification tests ==="

use "examples/texas_education.dta", clear

quietly nnsreg staar_pass spend_pupil, order(3) method(connect) partition(xy) ///
    noplot generate(chk_fit)
assert e(cmd) == "nnsreg"
assert e(N) == 250
assert e(order) == 3
assert "`e(method)'" == "connect"
assert "`e(model)'" == "connect"
assert "`e(vce)'" == "heuristic"
assert "`e(vcetype)'" == "Heuristic reference"
assert "`e(indepvar)'" == "spend_pupil"
assert "`e(marginsok)'" == ""
assert "`e(marginsnotok)'" == "_ALL"
assert "`e(estat_cmd)'" == "nnsreg_estat"
matrix chk_beta = e(beta)
local chk_beta_names : colnames chk_beta
assert "`: word 1 of `chk_beta_names''" == "spend_pupil"
quietly summarize staar_pass
local chk_ysd = sqrt(r(Var))
quietly summarize spend_pupil
local chk_xsd = sqrt(r(Var))
assert abs(chk_beta[1,1] - e(b)[1,1] * `chk_xsd' / `chk_ysd') < 1e-10
capture estat summarize
assert _rc == 321
assert "`e(knots)'" != ""
assert strpos("`e(cmdline)'", "nnsreg ") == 1

// predict must reproduce generate() exactly, and residuals must add up.
quietly predict double chk_p
quietly gen double chk_pdiff = abs(chk_p - chk_fit)
quietly summarize chk_pdiff
assert r(max) < 1e-8
quietly predict double chk_r, residuals
quietly gen double chk_rdiff = abs((staar_pass - chk_p) - chk_r)
quietly summarize chk_rdiff
assert r(max) < 1e-8
quietly gen double chk_r2 = chk_r^2
quietly summarize chk_r2, meanonly
assert abs(e(rss) - r(sum)) < 1e-8
assert missing(e(F)) & missing(e(ll)) & missing(e(ll_0))

// if/in restrictions respect the estimation sample.
quietly nnsreg staar_pass spend_pupil if spend_pupil < 15, order(2) noplot
assert e(N) < 250 & e(N) > 0

// Documented option abbreviations parse (ord, noi, meth, part, nop, gen).
quietly nnsreg staar_pass spend_pupil, ord(2) meth(ols) part(xonly) ///
    noi(median) nop gen(chk_abbrev)
assert e(order) == 2 & "`e(method)'" == "ols" & "`e(partition)'" == "xonly"

// Invalid options exit with error 198.
capture nnsreg staar_pass spend_pupil, method(bogus) noplot
assert _rc == 198
capture nnsreg staar_pass spend_pupil, noise(bogus) noplot
assert _rc == 198
capture nnsreg staar_pass spend_pupil, partition(bogus) noplot
assert _rc == 198
capture nnsreg staar_pass spend_pupil, order(0) noplot
assert _rc == 198

// vce() is deliberately restricted to robust and cluster for method(ols).
foreach bad in hc2 hc3 bootstrap jackknife {
    capture nnsreg staar_pass spend_pupil, method(ols) vce(`bad') noplot
    assert _rc == 198
}
capture nnsreg staar_pass spend_pupil, method(connect) vce(robust) noplot
assert _rc == 198
quietly nnsreg staar_pass spend_pupil, method(ols) vce(robust) noplot
assert "`e(vce)'" == "robust"
assert "`e(model)'" == "ols"
quietly nnsreg staar_pass spend_pupil, method(ols) vce(r) noplot
assert "`e(vce)'" == "robust"
capture margins, at(spend_pupil=(5 10 15))
assert _rc == 322

// Missing cluster IDs must be excluded before partitioning so all reported
// results describe the same estimation sample.
gen long chk_cluster = ceil(_n / 10)
replace chk_cluster = . in 1/10
quietly nnsreg staar_pass spend_pupil, method(ols) vce(cluster chk_cluster) ///
    slopes noplot generate(chk_cluster_fit)
assert e(N) == 240
assert e(N_clust) == 24
assert "`e(clustvar)'" == "chk_cluster"
count if e(sample)
assert r(N) == e(N)
count if !missing(chk_cluster_fit)
assert r(N) == e(N)
matrix chk_segments = e(segments)
mata: st_numscalar("chk_seg_n", sum(st_matrix("chk_segments")[.,5]))
assert chk_seg_n == e(N)

// order(1) is the unsplit full-sample partition.
quietly nnsreg staar_pass spend_pupil, order(1) method(ols) noplot
assert e(order) == 1
local order1_knots : word count `e(knots)'
quietly nnsreg staar_pass spend_pupil, order(2) method(ols) noplot
local order2_knots : word count `e(knots)'
assert `order2_knots' > `order1_knots'

di "Section 1 verification tests passed."

// ---------------------------------------------------------------------------
// 2. Education example: smooth S-curve with a plateau
//    Expectation (honest framing): the kernel smoother should recover a
//    smooth curve at least as well as NNS here.  The NNS payoff is the
//    stored slope components and the testable plateau (lincom below).
// ---------------------------------------------------------------------------
di _n "=== Section 2: education example (smooth S-curve) ==="

use "examples/texas_education.dta", clear

regress staar_pass spend_pupil
estimates store edu_linear
predict double fit_edu_linear, xb
calc_metrics, fitvar(fit_edu_linear) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("Linear regression") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_lin : di %5.3f r(rmse_true)

quietly npregress kernel staar_pass spend_pupil
predict double fit_edu_np, mean
calc_metrics, fitvar(fit_edu_np) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("npregress kernel") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_np : di %5.3f r(rmse_true)

quietly npregress series staar_pass spend_pupil
predict double fit_edu_series, mean
calc_metrics, fitvar(fit_edu_series) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("npregress series") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_series : di %5.3f r(rmse_true)

nnsreg staar_pass spend_pupil, order(2) method(connect) partition(xy) ///
    noplot generate(fit_edu_conn_o2)
estimates store edu_conn_o2
calc_metrics, fitvar(fit_edu_conn_o2) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("NNS connect xy o2") (r(r2)) (r(rmse_obs)) (r(rmse_true))

nnsreg staar_pass spend_pupil, order(3) method(connect) partition(xy) ///
    noplot generate(fit_edu_conn_o3)
estimates store edu_conn_o3
calc_metrics, fitvar(fit_edu_conn_o3) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("NNS connect xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_conn : di %5.3f r(rmse_true)

nnsreg staar_pass spend_pupil, order(3) method(connect) partition(xy) ///
    noise(median) noplot generate(fit_edu_conn_med_o3)
estimates store edu_conn_med_o3
calc_metrics, fitvar(fit_edu_conn_med_o3) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("NNS median xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))

nnsreg staar_pass spend_pupil, order(3) method(ols) partition(xy) ///
    noplot generate(fit_edu_ols_o3)
estimates store edu_ols_o3
calc_metrics, fitvar(fit_edu_ols_o3) yvar(staar_pass) truevar(true_staar)
post `diag' ("education") ("NNS OLS xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_ols : di %5.3f r(rmse_true)

// Figure: 2x2 small multiples, one fitted method per panel.
fitpanel, yvar(staar_pass) xvar(spend_pupil) truevar(true_staar) ///
    fitvar(fit_edu_linear) ptitle("A. Linear regression") ///
    psub("True-curve RMSE = `rt_lin'") gname(edu_p1) ylab(20(20)100)
fitpanel, yvar(staar_pass) xvar(spend_pupil) truevar(true_staar) ///
    fitvar(fit_edu_np) ptitle("B. npregress kernel") ///
    psub("True-curve RMSE = `rt_np'") gname(edu_p2) ylab(20(20)100)
fitpanel, yvar(staar_pass) xvar(spend_pupil) truevar(true_staar) ///
    fitvar(fit_edu_conn_o3) ptitle("C. NNS connect, order(3)") ///
    psub("True-curve RMSE = `rt_conn'") gname(edu_p3) ylab(20(20)100)
fitpanel, yvar(staar_pass) xvar(spend_pupil) truevar(true_staar) ///
    fitvar(fit_edu_ols_o3) ptitle("D. NNS OLS, order(3)") ///
    psub("True-curve RMSE = `rt_ols'") gname(edu_p4) ylab(20(20)100)

graph combine edu_p1 edu_p2 edu_p3 edu_p4, cols(2) imargin(small) ///
    l1title("STAAR Math Passing Rate (%)", size(small)) ///
    b1title("Per-Pupil Operational Spending (Thousands of Dollars)", size(small)) ///
    note("Gray dots: simulated observations.  Black line: true curve.  Orange line: fitted values.", ///
        size(vsmall)) ///
    graphregion(color(white)) xsize(7.5) ysize(6) name(edu_fits, replace)
graph export "manuscript/fig_edu_fits.pdf", replace
graph export "manuscript/fig_edu_fits.png", replace width(2200)

// Figure: stored slope components for the two NNS fits, with the adaptive
// knot locations written into the coefficient labels.  Restricting the plot
// to the two NNS models (dropping the linear fit) and renaming the
// coefficients removes the long-variable-label collisions of earlier drafts.
estimates restore edu_ols_o3
local clab `"spend_pupil = "Slope of first segment""'
local j = 1
foreach k in `e(knots)' {
    local kf = strtrim(string(`k', "%4.1f"))
    local clab `"`clab' knot`j' = "Slope change at spending = `kf'""'
    local ++j
}
coefplot ///
    (edu_conn_o3, label("NNS connect") mcolor("212 69 0") ciopts(lcolor("212 69 0"))) ///
    (edu_ols_o3, label("NNS OLS") mcolor("27 45 85") ciopts(lcolor("27 45 85"))), ///
    drop(_cons) coeflabels(`clab', labsize(small)) ///
    xline(0, lcolor(gs10) lpattern(dash)) ///
    xtitle("Slope or slope change (percentage points per thousand dollars)", ///
        size(small)) ///
    xlabel(, labsize(small)) ///
    legend(rows(1) size(small) position(6) region(lcolor(none))) ///
    graphregion(color(white)) xsize(6.5) ysize(5) name(edu_coef, replace)
graph export "manuscript/fig_edu_coef.pdf", replace
graph export "manuscript/fig_edu_coef.png", replace width(2000)

// Postestimation demo: standard Stata reporting tools work after nnsreg.
estimates table edu_conn_o3 edu_ols_o3, b(%9.3f) se(%9.3f) stats(r2 rmse N)

// Plateau test: the slope of the final segment equals the first-segment
// coefficient plus every knot (slope-change) coefficient.  A terminal slope
// near zero is the "plateau" read of the S-curve.  These standard errors are
// conditional on the NNS-selected knots, so read this as a descriptive check
// rather than exact inference.
estimates restore edu_ols_o3
local terms "spend_pupil"
foreach nm in `: colnames e(b)' {
    if strpos("`nm'", "knot") local terms "`terms' + `nm'"
}
di as text "Terminal-segment slope (plateau test), NNS OLS order(3):"
lincom `terms'
local tcrit = invttail(r(df), 0.025)

tempname fh
file open `fh' using "manuscript/table_postest.csv", write replace
file write `fh' "quantity,estimate,se,ci_lower,ci_upper" _n
file write `fh' "edu_terminal_slope,`r(estimate)',`r(se)'," ///
    "`=r(estimate) - `tcrit'*r(se)',`=r(estimate) + `tcrit'*r(se)'" _n

// Sensitivity figure: order and center choices, NNS connect variants only.
twoway ///
    (scatter staar_pass spend_pupil, mcolor("108 122 141%35") msymbol(Oh) msize(tiny)) ///
    (line true_staar spend_pupil, sort lcolor(black) lwidth(medthick)) ///
    (line fit_edu_conn_o2 spend_pupil, sort lcolor("43 108 176") lpattern(shortdash) lwidth(medthin)) ///
    (line fit_edu_conn_o3 spend_pupil, sort lcolor("212 69 0") lpattern(dash) lwidth(medthick)) ///
    (line fit_edu_conn_med_o3 spend_pupil, sort lcolor("0 109 44") lpattern(longdash) lwidth(medthin)), ///
    ytitle("STAAR Math Passing Rate (%)", size(small)) ///
    xtitle("Per-Pupil Operational Spending (Thousands of Dollars)", size(small)) ///
    ylabel(20(20)100, labsize(small) angle(horizontal)) xlabel(, labsize(small)) ///
    legend(order(2 "True curve" 3 "order(2), mean" 4 "order(3), mean" ///
        5 "order(3), median") rows(1) size(small) position(6) region(lcolor(none))) ///
    graphregion(color(white)) xsize(7) ysize(4.5) name(edu_sensitivity, replace)
graph export "manuscript/fig_edu_sensitivity.pdf", replace
graph export "manuscript/fig_edu_sensitivity.png", replace width(2000)

// ---------------------------------------------------------------------------
// 3. Water example: threshold behavior
//    The payoff here is interpretability: nnsreg stores its knots in
//    e(knots), so the analyst can display where the fitted slope changes and
//    compare those locations with substantive expectations.
// ---------------------------------------------------------------------------
di _n "=== Section 3: water example (threshold) ==="

use "examples/texas_water.dta", clear

regress water_loss pipe_age
estimates store water_linear
predict double fit_water_linear, xb
calc_metrics, fitvar(fit_water_linear) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("Linear regression") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_lin : di %5.3f r(rmse_true)

quietly npregress kernel water_loss pipe_age
predict double fit_water_np, mean
calc_metrics, fitvar(fit_water_np) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("npregress kernel") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_np : di %5.3f r(rmse_true)

quietly npregress series water_loss pipe_age
predict double fit_water_series, mean
calc_metrics, fitvar(fit_water_series) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("npregress series") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_series : di %5.3f r(rmse_true)

nnsreg water_loss pipe_age, order(2) method(connect) partition(xy) ///
    noplot generate(fit_water_conn_o2)
estimates store water_conn_o2
calc_metrics, fitvar(fit_water_conn_o2) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("NNS connect xy o2") (r(r2)) (r(rmse_obs)) (r(rmse_true))

nnsreg water_loss pipe_age, order(3) method(connect) partition(xy) ///
    noplot generate(fit_water_conn_o3)
estimates store water_conn_o3
calc_metrics, fitvar(fit_water_conn_o3) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("NNS connect xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_conn : di %5.3f r(rmse_true)

nnsreg water_loss pipe_age, order(3) method(connect) partition(xy) ///
    noise(median) noplot generate(fit_water_conn_med_o3)
estimates store water_conn_med_o3
calc_metrics, fitvar(fit_water_conn_med_o3) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("NNS median xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))

nnsreg water_loss pipe_age, order(3) method(ols) partition(xy) ///
    noplot generate(fit_water_ols_o3)
estimates store water_ols_o3
calc_metrics, fitvar(fit_water_ols_o3) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("NNS OLS xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_ols : di %5.3f r(rmse_true)

// partition(xonly) sensitivity row: vertical bands instead of quadrants.
nnsreg water_loss pipe_age, order(3) method(ols) partition(xonly) ///
    noplot generate(fit_water_xonly_o3)
estimates store water_xonly_o3
calc_metrics, fitvar(fit_water_xonly_o3) yvar(water_loss) truevar(true_loss)
post `diag' ("water") ("NNS OLS xonly o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))

// Figure: 2x2 small multiples, one fitted method per panel.
fitpanel, yvar(water_loss) xvar(pipe_age) truevar(true_loss) ///
    fitvar(fit_water_linear) ptitle("A. Linear regression") ///
    psub("True-curve RMSE = `rt_lin'") gname(water_p1) ylab(0(10)50)
fitpanel, yvar(water_loss) xvar(pipe_age) truevar(true_loss) ///
    fitvar(fit_water_np) ptitle("B. npregress kernel") ///
    psub("True-curve RMSE = `rt_np'") gname(water_p2) ylab(0(10)50)
fitpanel, yvar(water_loss) xvar(pipe_age) truevar(true_loss) ///
    fitvar(fit_water_conn_o3) ptitle("C. NNS connect, order(3)") ///
    psub("True-curve RMSE = `rt_conn'") gname(water_p3) ylab(0(10)50)
fitpanel, yvar(water_loss) xvar(pipe_age) truevar(true_loss) ///
    fitvar(fit_water_ols_o3) ptitle("D. NNS OLS, order(3)") ///
    psub("True-curve RMSE = `rt_ols'") gname(water_p4) ylab(0(10)50)

graph combine water_p1 water_p2 water_p3 water_p4, cols(2) imargin(small) ///
    l1title("Annual Water Loss Rate (%)", size(small)) ///
    b1title("Average Water Pipe Age (Years)", size(small)) ///
    note("Gray dots: simulated observations.  Black line: true curve.  Orange line: fitted values.", ///
        size(vsmall)) ///
    graphregion(color(white)) xsize(7.5) ysize(6) name(water_fits, replace)
graph export "manuscript/fig_water_fits.pdf", replace
graph export "manuscript/fig_water_fits.png", replace width(2200)

// Figure: where did nnsreg place its knots?  The e(knots) macro makes the
// slope-change locations reportable; here they are drawn against the true
// change point at age 35 that was built into the simulation.
estimates restore water_ols_o3
local klist ""
foreach k in `e(knots)' {
    local klist "`klist' `k'"
}
twoway ///
    (scatter water_loss pipe_age, mcolor("108 122 141%35") msymbol(Oh) msize(tiny)) ///
    (line true_loss pipe_age, sort lcolor(black) lwidth(medthin)) ///
    (line fit_water_ols_o3 pipe_age, sort lcolor("212 69 0") lwidth(medthick)), ///
    xline(`klist', lcolor("108 122 141") lpattern(shortdash) lwidth(vthin)) ///
    xline(35, lcolor("27 45 85") lwidth(medthick)) ///
    text(46 34 "true change point (age 35)", placement(w) box ///
        fcolor(white) lcolor(none) margin(small) size(small) ///
        color("27 45 85")) ///
    ytitle("Annual Water Loss Rate (%)", size(small)) ///
    xtitle("Average Water Pipe Age (Years)", size(small)) ///
    ylabel(0(10)50, labsize(small) angle(horizontal)) xlabel(, labsize(small)) ///
    legend(order(2 "True curve" 3 "NNS OLS fit, order(3)") rows(1) ///
        size(small) position(6) region(lcolor(none))) ///
    note("Vertical dashed gray lines: NNS interior knots from e(knots).  Navy vertical line: true change point.", ///
        size(vsmall)) ///
    graphregion(color(white)) xsize(7) ysize(4.5) name(water_knots, replace)
graph export "manuscript/fig_water_knots.pdf", replace
graph export "manuscript/fig_water_knots.png", replace width(2000)

// Slope contrast: first segment (young pipes) vs. terminal segment (old
// pipes).  Same conditional-on-knots caveat as the education plateau test.
di as text "First-segment slope (young pipes), NNS OLS order(3):"
lincom pipe_age
local tcrit = invttail(r(df), 0.025)
file write `fh' "water_first_slope,`r(estimate)',`r(se)'," ///
    "`=r(estimate) - `tcrit'*r(se)',`=r(estimate) + `tcrit'*r(se)'" _n

local terms "pipe_age"
foreach nm in `: colnames e(b)' {
    if strpos("`nm'", "knot") local terms "`terms' + `nm'"
}
di as text "Terminal-segment slope (old pipes), NNS OLS order(3):"
lincom `terms'
local tcrit = invttail(r(df), 0.025)
file write `fh' "water_terminal_slope,`r(estimate)',`r(se)'," ///
    "`=r(estimate) - `tcrit'*r(se)',`=r(estimate) + `tcrit'*r(se)'" _n
file close `fh'

// Sensitivity figure: order and center choices, NNS connect variants only.
twoway ///
    (scatter water_loss pipe_age, mcolor("108 122 141%35") msymbol(Oh) msize(tiny)) ///
    (line true_loss pipe_age, sort lcolor(black) lwidth(medthick)) ///
    (line fit_water_conn_o2 pipe_age, sort lcolor("43 108 176") lpattern(shortdash) lwidth(medthin)) ///
    (line fit_water_conn_o3 pipe_age, sort lcolor("212 69 0") lpattern(dash) lwidth(medthick)) ///
    (line fit_water_conn_med_o3 pipe_age, sort lcolor("0 109 44") lpattern(longdash) lwidth(medthin)), ///
    ytitle("Annual Water Loss Rate (%)", size(small)) ///
    xtitle("Average Water Pipe Age (Years)", size(small)) ///
    ylabel(0(10)50, labsize(small) angle(horizontal)) xlabel(, labsize(small)) ///
    legend(order(2 "True curve" 3 "order(2), mean" 4 "order(3), mean" ///
        5 "order(3), median") rows(1) size(small) position(6) region(lcolor(none))) ///
    graphregion(color(white)) xsize(7) ysize(4.5) name(water_sensitivity, replace)
graph export "manuscript/fig_water_sensitivity.pdf", replace
graph export "manuscript/fig_water_sensitivity.png", replace width(2000)

// ---------------------------------------------------------------------------
// 4. Clustered peak/dip benchmark: the case for segmentation
//    Data are dense in three regions with a sharp dip between them, the
//    setting where Vinod and Viole (2018) report partition-based fits doing
//    better than a single global bandwidth.
// ---------------------------------------------------------------------------
di _n "=== Section 4: clustered peak/dip benchmark ==="

use "examples/clustered_peak_dip.dta", clear

regress outcome_index cluster_score
estimates store cluster_linear
predict double fit_cluster_linear, xb
calc_metrics, fitvar(fit_cluster_linear) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("Linear regression") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_lin : di %5.3f r(rmse_true)

quietly npregress kernel outcome_index cluster_score
predict double fit_cluster_np, mean
calc_metrics, fitvar(fit_cluster_np) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("npregress kernel") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_np : di %5.3f r(rmse_true)

quietly npregress series outcome_index cluster_score
predict double fit_cluster_series, mean
calc_metrics, fitvar(fit_cluster_series) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("npregress series") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_series : di %5.3f r(rmse_true)

nnsreg outcome_index cluster_score, order(3) method(connect) partition(xy) ///
    noplot generate(fit_cluster_conn_o3)
estimates store cluster_conn_o3
calc_metrics, fitvar(fit_cluster_conn_o3) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("NNS connect xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_conn : di %5.3f r(rmse_true)

nnsreg outcome_index cluster_score, order(3) method(ols) partition(xy) ///
    noplot generate(fit_cluster_ols_o3)
estimates store cluster_ols_o3
calc_metrics, fitvar(fit_cluster_ols_o3) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("NNS OLS xy o3") (r(r2)) (r(rmse_obs)) (r(rmse_true))
local rt_ols : di %5.3f r(rmse_true)

nnsreg outcome_index cluster_score, order(4) method(ols) partition(xy) ///
    noplot generate(fit_cluster_ols_o4)
estimates store cluster_ols_o4
calc_metrics, fitvar(fit_cluster_ols_o4) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("NNS OLS xy o4") (r(r2)) (r(rmse_obs)) (r(rmse_true))

nnsreg outcome_index cluster_score, order(5) method(ols) partition(xy) ///
    noplot generate(fit_cluster_ols_o5)
estimates store cluster_ols_o5
calc_metrics, fitvar(fit_cluster_ols_o5) yvar(outcome_index) truevar(true_index)
post `diag' ("clustered") ("NNS OLS xy o5") (r(r2)) (r(rmse_obs)) (r(rmse_true))

// Figure: 2x2 small multiples, one fitted method per panel.
fitpanel, yvar(outcome_index) xvar(cluster_score) truevar(true_index) ///
    fitvar(fit_cluster_linear) ptitle("A. Linear regression") ///
    psub("True-curve RMSE = `rt_lin'") gname(cl_p1) ylab(0(10)40)
fitpanel, yvar(outcome_index) xvar(cluster_score) truevar(true_index) ///
    fitvar(fit_cluster_np) ptitle("B. npregress kernel") ///
    psub("True-curve RMSE = `rt_np'") gname(cl_p2) ylab(0(10)40)
fitpanel, yvar(outcome_index) xvar(cluster_score) truevar(true_index) ///
    fitvar(fit_cluster_conn_o3) ptitle("C. NNS connect, order(3)") ///
    psub("True-curve RMSE = `rt_conn'") gname(cl_p3) ylab(0(10)40)
fitpanel, yvar(outcome_index) xvar(cluster_score) truevar(true_index) ///
    fitvar(fit_cluster_ols_o3) ptitle("D. NNS OLS, order(3)") ///
    psub("True-curve RMSE = `rt_ols'") gname(cl_p4) ylab(0(10)40)

graph combine cl_p1 cl_p2 cl_p3 cl_p4, cols(2) imargin(small) ///
    l1title("Outcome Index", size(small)) ///
    b1title("Clustered Index Score", size(small)) ///
    note("Gray dots: simulated observations.  Black line: true curve.  Orange line: fitted values.", ///
        size(vsmall)) ///
    graphregion(color(white)) xsize(7.5) ysize(6) name(cluster_fits, replace)
graph export "manuscript/fig_cluster_fits.pdf", replace
graph export "manuscript/fig_cluster_fits.png", replace width(2200)

// Sensitivity figure: partition depth (order) against the kernel benchmark.
twoway ///
    (scatter outcome_index cluster_score, mcolor("108 122 141%35") msymbol(Oh) msize(tiny)) ///
    (line true_index cluster_score, sort lcolor(black) lwidth(medthick)) ///
    (line fit_cluster_np cluster_score, sort lcolor("108 122 141") lpattern(dot) lwidth(medthin)) ///
    (line fit_cluster_ols_o3 cluster_score, sort lcolor("212 69 0") lpattern(dash) lwidth(medthick)) ///
    (line fit_cluster_ols_o4 cluster_score, sort lcolor("43 108 176") lpattern(longdash) lwidth(medthin)) ///
    (line fit_cluster_ols_o5 cluster_score, sort lcolor("0 109 44") lpattern(shortdash) lwidth(medthin)), ///
    ytitle("Outcome Index", size(small)) ///
    xtitle("Clustered Index Score", size(small)) ///
    ylabel(0(10)40, labsize(small) angle(horizontal)) xlabel(, labsize(small)) ///
    legend(order(2 "True curve" 3 "npregress" 4 "NNS OLS(3)" 5 "NNS OLS(4)" ///
        6 "NNS OLS(5)") rows(1) size(small) position(6) region(lcolor(none))) ///
    graphregion(color(white)) xsize(7) ysize(4.5) name(cluster_sensitivity, replace)
graph export "manuscript/fig_cluster_sensitivity.pdf", replace
graph export "manuscript/fig_cluster_sensitivity.png", replace width(2000)

// ---------------------------------------------------------------------------
// 5. Monte Carlo study: repeat the comparison across `MCREPS' fresh draws
//    A single simulated dataset can flatter any method by luck of the draw.
//    Here every design is redrawn `MCREPS' times (new x values and new
//    noise), each method is refit, and RMSE against the known true curve is
//    recorded.  Seeds are 10000 + 3*(rep-1) + design index, so every draw is
//    reproducible.
// ---------------------------------------------------------------------------
di _n "=== Section 5: Monte Carlo study (`MCREPS' replications) ==="

tempname mc
postfile `mc' str14 design int rep str8 mcode str24 model ///
    double(r2 rmse_obs rmse_true) using "examples/mc_results.dta", replace

local np_failures = 0
local series_failures = 0
local nns_failures = 0
timer clear 1
timer on 1
forvalues r = 1/`MCREPS' {
    local dnum = 0
    foreach d in education water clustered {
        local ++dnum
        if "`d'" == "education" {
            local y staar_pass
            local x spend_pupil
            local t true_staar
        }
        else if "`d'" == "water" {
            local y water_loss
            local x pipe_age
            local t true_loss
        }
        else {
            local y outcome_index
            local x cluster_score
            local t true_index
        }

        nns_simdata, design(`d') seed(`=10000 + 3*(`r'-1) + `dnum'')

        quietly regress `y' `x'
        quietly predict double f_lin, xb
        calc_metrics, fitvar(f_lin) yvar(`y') truevar(`t')
        post `mc' ("`d'") (`r') ("linear") ("Linear regression") ///
            (r(r2)) (r(rmse_obs)) (r(rmse_true))

        capture {
            quietly npregress kernel `y' `x'
            quietly predict double f_np, mean
        }
        if _rc == 0 {
            calc_metrics, fitvar(f_np) yvar(`y') truevar(`t')
            post `mc' ("`d'") (`r') ("kernel") ("npregress kernel") ///
                (r(r2)) (r(rmse_obs)) (r(rmse_true))
        }
        else {
            local ++np_failures
            post `mc' ("`d'") (`r') ("kernel") ("npregress kernel") (.) (.) (.)
        }

        capture {
            quietly npregress series `y' `x'
            quietly predict double f_series, mean
        }
        if _rc == 0 {
            calc_metrics, fitvar(f_series) yvar(`y') truevar(`t')
            post `mc' ("`d'") (`r') ("series") ("npregress series") ///
                (r(r2)) (r(rmse_obs)) (r(rmse_true))
        }
        else {
            local ++series_failures
            post `mc' ("`d'") (`r') ("series") ("npregress series") (.) (.) (.)
        }

        local i = 0
        foreach spec in "order(3) method(connect)" "order(3) method(ols)" ///
            "order(4) method(ols)" {
            local ++i
            local mcode : word `i' of conn3 ols3 ols4
            local mname : word `i' of "NNS connect, order(3)" ///
                "NNS OLS, order(3)" "NNS OLS, order(4)"
            capture {
                quietly nnsreg `y' `x', `spec' partition(xy) ///
                    noplot generate(f_nns`i')
            }
            if _rc == 0 {
                calc_metrics, fitvar(f_nns`i') yvar(`y') truevar(`t')
                post `mc' ("`d'") (`r') ("`mcode'") (`"`mname'"') ///
                    (r(r2)) (r(rmse_obs)) (r(rmse_true))
            }
            else {
                local ++nns_failures
                post `mc' ("`d'") (`r') ("`mcode'") (`"`mname'"') (.) (.) (.)
                di as error "  nnsreg failed: rep `r', design `d', `spec'"
            }
        }
    }
    if mod(`r', 10) == 0 di as text "  Monte Carlo replication `r' of `MCREPS'"
}
timer off 1
postclose `mc'
timer list 1
di as text "npregress kernel failures across all replications: `np_failures'"
di as text "npregress series failures across all replications: `series_failures'"
di as text "nnsreg failures across all replications: `nns_failures'"
// Certification: nnsreg must fit every replication draw without error.
assert `nns_failures' == 0
assert `np_failures' == 0
assert `series_failures' == 0

// Summary table: mean and SD of true-curve RMSE by design and model, plus
// the share of replications in which each model beats npregress kernel.
use "examples/mc_results.dta", clear
preserve
keep if mcode == "kernel"
keep design rep rmse_true
rename rmse_true rmse_true_kernel
tempfile kern
save `kern'
restore
merge m:1 design rep using `kern', nogenerate
gen byte beats_kernel = rmse_true < rmse_true_kernel ///
    if !missing(rmse_true, rmse_true_kernel) & mcode != "kernel"

gen msort = cond(mcode == "linear", 1, cond(mcode == "kernel", 2, ///
    cond(mcode == "series", 3, cond(mcode == "conn3", 4, ///
    cond(mcode == "ols3", 5, 6)))))

collapse (mean) mean_rmse_true = rmse_true (sd) sd_rmse_true = rmse_true ///
    (mean) mean_rmse_obs = rmse_obs (mean) share_beats_kernel = beats_kernel ///
    (count) reps = rmse_true, by(design mcode model msort)
sort design msort
format mean_rmse_true sd_rmse_true mean_rmse_obs %6.3f
format share_beats_kernel %5.2f
list design model mean_rmse_true sd_rmse_true share_beats_kernel reps, ///
    sepby(design) noobs
export delimited design model mean_rmse_true sd_rmse_true mean_rmse_obs ///
    share_beats_kernel reps using "manuscript/table_mc_summary.csv", replace

// Figure: distribution of true-curve RMSE across replications.  Linear
// regression is omitted (its RMSE is several times larger and would
// compress the boxes); it appears in the summary table instead.
use "examples/mc_results.dta", clear
drop if mcode == "linear" | missing(rmse_true)
gen mnum = cond(mcode == "kernel", 1, cond(mcode == "series", 2, ///
    cond(mcode == "conn3", 3, cond(mcode == "ols3", 4, 5))))
label define mnum 1 "np kernel" 2 "np series" 3 "NNS conn(3)" ///
    4 "NNS OLS(3)" 5 "NNS OLS(4)"
label values mnum mnum
gen dnum = cond(design == "education", 1, cond(design == "water", 2, 3))
label define dnum 1 "Education: smooth S-curve" 2 "Water: threshold" ///
    3 "Clustered: peak and dip"
label values dnum dnum

graph box rmse_true, over(mnum, label(angle(45) labsize(small))) ///
    by(dnum, yrescale cols(3) note("") graphregion(color(white))) ///
    ytitle("RMSE against the true curve", size(small)) ///
    ylabel(, labsize(small) angle(horizontal)) ///
    box(1, color("27 45 85")) ///
    xsize(9) ysize(4) name(mc_rmse, replace)
graph export "manuscript/fig_mc_rmse.pdf", replace
graph export "manuscript/fig_mc_rmse.png", replace width(2400)

// ---------------------------------------------------------------------------
// 6. Single-run diagnostics table (manuscript table 1)
// ---------------------------------------------------------------------------
postclose `diag'

use "examples/model_diagnostics.dta", clear
format r2 %6.3f
format rmse_obs rmse_true %6.3f
export delimited using "manuscript/table_model_diagnostics.csv", replace
list, sepby(example)

di _n "All tests, comparisons, sensitivity checks, and figures completed successfully."
