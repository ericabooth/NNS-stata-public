*! version 1.0.0  18jul2026
*! Two robustness checks for the nnsreg Stata Journal article, kept separate
*! from the main replication script so that its outputs are unchanged:
*!
*!   1. Empirical coverage of the nominal-95% conditional interval for the
*!      terminal-segment slope, across 100 replications per simulated design,
*!      for method(ols) and method(connect).  The target is the true secant
*!      slope of the generating curve over the fitted terminal segment.  Poor
*!      coverage bounds the inferential use of these intervals and supports the
*!      "descriptive, not tests" reading.
*!
*!   2. An alternative clustered design (noise SD raised from 1.2 to 2.0) to
*!      test whether the NNS OLS between-replication variance advantage
*!      persists, and a percentile-knot linear spline (mkspline) arm to show
*!      the adaptive-knot advantage on data with known truth.
*!
*! Writes manuscript/table_coverage.csv and manuscript/table_altclust.csv.
clear all
set more off
adopath + "code"

// Defines nns_simdata (also rewrites the canonical datasets; harmless here).
do "examples/simulated_data.do"

// True generating curve evaluated at an arbitrary x, by design.
capture program drop nns_true_at
program define nns_true_at, rclass
    syntax, DESign(string) X(real)
    if "`design'" == "education" {
        return scalar y = 40 + 45 / (1 + exp(-0.8 * (`x' - 11)))
    }
    else if "`design'" == "water" {
        if `x' <= 35 return scalar y = 8
        else return scalar y = 8 + 0.4 * (`x' - 35) + 0.01 * (`x' - 35)^2
    }
    else {
        return scalar y = 18 + 0.18 * `x' ///
            + 10 * exp(-((`x' - 24)^2) / 35) ///
            - 16 * exp(-((`x' - 35)^2) / 18) ///
            + 8 * exp(-((`x' - 58)^2) / 50)
    }
end

// ---------------------------------------------------------------------------
// 1. Coverage of the conditional terminal-segment interval
// ---------------------------------------------------------------------------
di _n "=== Coverage check: 100 replications per design ==="
tempname cov
postfile `cov' str12 design str8 method byte covered ///
    using "examples/coverage_raw.dta", replace

forvalues r = 1/100 {
    local dnum = 0
    foreach d in education water clustered {
        local ++dnum
        if "`d'" == "education" {
            local y staar_pass
            local x spend_pupil
        }
        else if "`d'" == "water" {
            local y water_loss
            local x pipe_age
        }
        else {
            local y outcome_index
            local x cluster_score
        }
        nns_simdata, design(`d') seed(`=20000 + 3*(`r'-1) + `dnum'')
        foreach m in ols connect {
            capture {
                quietly nnsreg `y' `x', order(3) method(`m') noplot
                matrix S = e(segments)
                local ns = rowsof(S)
                local xlo = S[`ns', 1]
                local xhi = S[`ns', 2]
                local slope = S[`ns', 3]
                local se = S[`ns', 4]
                local tc = invttail(e(df_r), 0.025)
                local lo = `slope' - `tc' * `se'
                local hi = `slope' + `tc' * `se'
                nns_true_at, design(`d') x(`xlo')
                local ylo = r(y)
                nns_true_at, design(`d') x(`xhi')
                local yhi = r(y)
                local truesec = (`yhi' - `ylo') / (`xhi' - `xlo')
                local cvd = (`lo' <= `truesec') & (`truesec' <= `hi')
                post `cov' ("`d'") ("`m'") (`cvd')
            }
        }
    }
    if mod(`r', 25) == 0 di as text "  coverage replication `r' of 100"
}
postclose `cov'

use "examples/coverage_raw.dta", clear
collapse (mean) coverage = covered (count) n = covered, by(design method)
gen dnum = cond(design == "education", 1, cond(design == "water", 2, 3))
sort dnum method
format coverage %5.3f
list design method coverage n, sepby(design) noobs
export delimited design method coverage n ///
    using "manuscript/table_coverage.csv", replace

// ---------------------------------------------------------------------------
// 2. Alternative clustered design (higher noise) + percentile-knot arm
// ---------------------------------------------------------------------------
di _n "=== Alternative clustered design: noise SD 2.0, 100 replications ==="
tempname alt
postfile `alt' str28 model double rmse_true int rep ///
    using "examples/altclust_raw.dta", replace

capture program drop nns_rmse_true
program define nns_rmse_true, rclass
    syntax, Fitvar(name) Truevar(name)
    tempvar e2
    quietly gen double `e2' = (`truevar' - `fitvar')^2
    quietly summarize `e2', meanonly
    return scalar rmse = sqrt(r(mean))
end

forvalues r = 1/100 {
    clear
    set seed `=30000 + `r''
    set obs 260
    gen mix = runiform()
    gen cluster_score = cond(mix < .45, rnormal(22, 3), ///
        cond(mix < .75, rnormal(36, 2.2), rnormal(58, 4)))
    replace cluster_score = max(5, min(70, cluster_score))
    gen true_index = 18 + .18 * cluster_score ///
        + 10 * exp(-((cluster_score - 24)^2) / 35) ///
        - 16 * exp(-((cluster_score - 35)^2) / 18) ///
        + 8 * exp(-((cluster_score - 58)^2) / 50)
    gen outcome_index = true_index + rnormal(0, 2.0)
    drop mix

    quietly npregress kernel outcome_index cluster_score
    quietly predict double f_k, mean
    nns_rmse_true, fitvar(f_k) truevar(true_index)
    post `alt' ("npregress kernel") (r(rmse)) (`r')
    drop f_k

    quietly npregress series outcome_index cluster_score
    quietly predict double f_s, mean
    nns_rmse_true, fitvar(f_s) truevar(true_index)
    post `alt' ("npregress series") (r(rmse)) (`r')
    drop f_s

    quietly nnsreg outcome_index cluster_score, order(3) method(ols) noplot
    quietly predict double f_o
    nns_rmse_true, fitvar(f_o) truevar(true_index)
    post `alt' ("NNS OLS, order(3)") (r(rmse)) (`r')
    drop f_o

    // Percentile-knot linear spline (rank-placed knots), 6 pieces.
    capture {
        mkspline sp 6 = cluster_score, pctile
        quietly regress outcome_index sp*
        quietly predict double f_m, xb
        nns_rmse_true, fitvar(f_m) truevar(true_index)
        post `alt' ("Percentile-knot spline") (r(rmse)) (`r')
        drop f_m sp*
    }
    if mod(`r', 25) == 0 di as text "  alt-clustered replication `r' of 100"
}
postclose `alt'

use "examples/altclust_raw.dta", clear
gen msort = cond(model == "npregress kernel", 1, ///
    cond(model == "npregress series", 2, ///
    cond(model == "NNS OLS, order(3)", 3, 4)))
collapse (mean) mean_rmse_true = rmse_true (sd) sd_rmse_true = rmse_true ///
    (count) reps = rmse_true, by(model msort)
sort msort
format mean_rmse_true sd_rmse_true %6.3f
list model mean_rmse_true sd_rmse_true reps, noobs
export delimited model mean_rmse_true sd_rmse_true reps ///
    using "manuscript/table_altclust.csv", replace

di _n "ROBUSTNESS CHECKS COMPLETE"
