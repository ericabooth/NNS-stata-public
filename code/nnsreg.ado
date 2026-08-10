*! version 1.3.2  18jul2026
*! nnsreg: Nonlinear Nonparametric Statistics using Partial Moments in Stata
*! Author: Eric A. Booth (Senior Researcher, Texas 2036, eric.a.booth@gmail.com)
*! 1.3.2: label connect-fit covariance/model metadata as heuristic/connect;
*!        remove inherited regress postestimation flags that nnsreg does not
*!        support; reject unsupported estat/margins calls cleanly; refresh
*!        standardized coefficients after reposting; trim the stored
*!        independent-variable name.
*! 1.3.1: restrict vce() to robust/cluster; align cluster-missingness across
*!        partitioning, estimation, prediction, and segment counts; clarify
*!        order-depth implementation; warn on empty and thin x ranges; align
*!        all retained connect-fit statistics with connected predictions.
*! 1.3.0: new options minobs(), level(), vce() (method(ols) only), slopes,
*!        knotlines; new stored results e(segments), e(knotmat), e(xmin),
*!        e(xmax), e(minobs); guards for degenerate (near-constant) regressors
*!        and perfect-fit covariance scaling; generate() now validated at parse
*!        time and confirmed before estimation; predict notes out-of-range
*!        scoring.
*! 1.2.1: guard against missing regression-point coordinates when the
*!        regressor support is clustered or gapped (empty central interval);
*!        method(connect) e(r2), e(rmse), and e(V) scaling now use raw
*!        (uncentered) residual sums of squares; option abbreviations
*!        (ord, noi, meth, part, nop, gen) now accepted as documented

program define nnsreg, eclass
    version 16.0
    local nns_cmdline "nnsreg `0'"
    syntax varlist(min=2 max=2) [if] [in] [, ///
        ORDer(integer 3) ///
        NOIse(string) ///
        METHod(string) ///
        PARTition(string) ///
        MINobs(integer 8) ///
        Level(cilevel) ///
        vce(passthru) ///
        SLOPes ///
        NOPlot ///
        KNOTLines ///
        graph_opts(string) ///
        GENerate(name) ]

    // 1. Validate options
    if "`noise'" == "" local noise "mean"
    if "`method'" == "" local method "connect"
    if "`partition'" == "" local partition "xy"

    if !inlist("`noise'", "mean", "median") {
        di as error "noise() must be either mean or median"
        exit 198
    }
    if !inlist("`method'", "connect", "ols") {
        di as error "method() must be either connect or ols"
        exit 198
    }
    if !inlist("`partition'", "xy", "xonly") {
        di as error "partition() must be either xy or xonly"
        exit 198
    }
    if `order' < 1 {
        di as error "order() must be a positive integer"
        exit 198
    }
    if `minobs' < 2 {
        di as error "minobs() must be at least 2"
        exit 198
    }
    // vce() requests a robust or cluster-robust covariance estimate and only
    // composes with method(ols); the connect fit posts a rescaled OLS
    // covariance that a sandwich estimator would not compose with, so refuse
    // the combination.
    if `"`vce'"' != "" & "`method'" == "connect" {
        di as error "vce() is allowed only with method(ols)"
        exit 198
    }

    // vce(passthru) deliberately leaves parsing to this command.  Restrict it
    // to the two documented choices so resampling VCEs cannot be mistaken for
    // valid inference with data-derived knots.
    local vcetype ""
    local clustvar ""
    if `"`vce'"' != "" {
        local vcebody = substr(`"`vce'"', 5, length(`"`vce'"') - 5)
        gettoken vcetype vcerest : vcebody
        local vcetype_l = lower("`vcetype'")
        local isrobust = strlen("`vcetype_l'") >= 1 & ///
            substr("robust", 1, strlen("`vcetype_l'")) == "`vcetype_l'"
        local iscluster = strlen("`vcetype_l'") >= 2 & ///
            substr("cluster", 1, strlen("`vcetype_l'")) == "`vcetype_l'"
        if `isrobust' {
            if `"`vcerest'"' != "" {
                di as error "vce(robust) does not take an additional argument"
                exit 198
            }
        }
        else if `iscluster' {
            gettoken clustvar vceextra : vcerest
            if "`clustvar'" == "" | `"`vceextra'"' != "" {
                di as error "vce(cluster) requires one cluster variable"
                exit 198
            }
            confirm variable `clustvar'
        }
        else {
            di as error "vce() must be robust or cluster clustvar"
            exit 198
        }
    }

    // 2. Mark the estimation sample
    marksample touse
    // The partition, spline regression, stored range, generated values, and
    // segment counts must use the same observations.  regress drops missing
    // cluster identifiers, so remove them before the partition is built.
    if "`clustvar'" != "" markout `touse' `clustvar'
    gettoken depvar indepvar : varlist
    local indepvar = strtrim("`indepvar'")

    // Confirm the generate() target now, before any estimation work, so a
    // name clash fails fast instead of after the fit.
    if "`generate'" != "" confirm new variable `generate'

    // Count observations
    quietly count if `touse'
    local N = r(N)
    if `N' == 0 {
        error 2000
    }
    if `N' < 3 {
        di as error "nnsreg needs at least 3 observations"
        exit 2001
    }

    // 3. Run partitioning and compute knots and NNS coefficients in Mata
    tempname b_nns knotmat
    local knots ""
    local nnsfail ""
    mata: nns_reg_mata("`depvar'", "`indepvar'", "`touse'", `order', `minobs', "`noise'", "`partition'", "`b_nns'", "knots", "`knotmat'", "nnsfail")

    // Degenerate support (near-constant regressor, or a partition that
    // collapses to a single point) is reported by the Mata core as a
    // friendly message rather than an opaque subscript error.
    if "`nnsfail'" != "" {
        di as error "`nnsfail'"
        exit 198
    }

    // 4. Generate spline basis variables in Stata
    local num_knots : word count `knots'

    // Temporary variables for spline basis
    tempvar sp_1
    quietly gen double `sp_1' = `indepvar' if `touse'
    local sp_vars `sp_1'
    local coef_names `indepvar'

    forvalues j = 1/`num_knots' {
        local k : word `j' of `knots'
        tempvar sp_`=`j'+1'
        quietly gen double `sp_`=`j'+1'' = cond(`indepvar' > `k', `indepvar' - `k', 0) if `touse'
        local sp_vars `sp_vars' `sp_`=`j'+1''
        local coef_names `coef_names' knot`j'
    }

    // 5. Fit spline model by OLS to set up e() structure.  For method(ols)
    // any user vce() is applied here; for method(connect) the OLS fit is only
    // scaffolding for the basis and degrees of freedom, so it stays plain.
    if "`method'" == "ols" {
        quietly regress `depvar' `sp_vars' if `touse', `vce'
    }
    else {
        quietly regress `depvar' `sp_vars' if `touse'
    }
    local df_m = e(df_m)
    local df_r = e(df_r)
    local r2_ols = e(r2)
    local rmse_ols = e(rmse)

    // A saturated spline leaves no residual degrees of freedom; stop with a
    // clear message rather than dividing by zero below.
    if `df_r' <= 0 {
        di as error "the spline saturates the data (no residual degrees of freedom); increase N, lower order(), or raise minobs()"
        exit 198
    }

    // 6. Handle estimation methods
    if "`method'" == "connect" {
        // Set column names for the NNS coefficient matrix
        matrix colnames `b_nns' = `coef_names' _cons

        // Calculate NNS fitted values and scale e(V)
        tempvar y_hat res
        quietly predict double `y_hat' if `touse'
        // Overwrite standard predictions temporarily with NNS connect values
        // to calculate NNS residuals
        mata: nns_calc_fitted("`sp_vars'", "`b_nns'", "`touse'", "`y_hat'")
        quietly gen double `res' = `depvar' - `y_hat' if `touse'

        // Use raw (uncentered) sums of squared residuals: the connected fit
        // is not least squares, so its residuals need not have mean zero and
        // a centered variance would overstate fit
        tempvar res2
        quietly gen double `res2' = `res'^2 if `touse'
        quietly sum `res2' if `touse', meanonly
        local rss = r(sum)
        local s2_nns = `rss' / `df_r'
        local s2_ols = `rmse_ols'^2

        // Scale the covariance matrix V.  A perfect OLS spline fit leaves
        // s2_ols == 0, which would make the scale factor missing; keep the
        // (degenerate, all-zero) covariance and warn instead of crashing.
        tempname V_nns
        if `s2_ols' == 0 {
            di as text "note: the linear-spline fit is perfect in sample; reported standard errors are degenerate"
            matrix `V_nns' = e(V)
        }
        else {
            matrix `V_nns' = e(V) * (`s2_nns' / `s2_ols')
        }

        // Repost the NNS coefficients and scaled V
        ereturn repost b = `b_nns' V = `V_nns', rename

        // Re-compute R-squared for NNS connect fit
        quietly sum `depvar' if `touse'
        local tss = r(Var) * (`N' - 1)
        if `tss' == 0 local r2 = .
        else local r2 = 1 - (`rss' / `tss')

        ereturn scalar r2 = `r2'
        ereturn scalar rmse = sqrt(`s2_nns')
        ereturn scalar rss = `rss'
        ereturn scalar mss = `tss' - `rss'
        ereturn scalar r2_a = 1 - (`rss' / `df_r') / (`tss' / (`N' - 1))
        // The OLS scaffold's F statistic and likelihood values do not
        // describe the connected coefficients.  Overwrite them rather than
        // leaving stale regress results available to collect/user code.
        ereturn scalar F = .
        ereturn scalar ll = .
        ereturn scalar ll_0 = .
        // Do not retain the OLS scaffold's method labels.  The posted
        // covariance is a heuristic reference for the connected fit.
        ereturn local vce "heuristic"
        ereturn local vcetype "Heuristic reference"
        ereturn local model "connect"
    }
    else {
        // method(ols) - keep the standard OLS coefficients and V
        local r2 = `r2_ols'
        tempname b_ols
        matrix `b_ols' = e(b)
        matrix colnames `b_ols' = `coef_names' _cons
        ereturn repost b = `b_ols', rename
    }

    // regress stores standardized coefficients in e(beta).  Reposting e(b)
    // does not refresh that matrix, so rebuild it for the final coefficients
    // and apply the public coefficient names.
    tempname beta_std
    local ncoef = `num_knots' + 1
    matrix `beta_std' = J(1, `ncoef', .)
    quietly summarize `depvar' if `touse'
    local ysd = sqrt(r(Var))
    forvalues j = 1/`ncoef' {
        local spj : word `j' of `sp_vars'
        quietly summarize `spj' if `touse'
        if `ysd' > 0 & r(Var) > 0 {
            matrix `beta_std'[1, `j'] = e(b)[1, `j'] * sqrt(r(Var)) / `ysd'
        }
    }
    matrix colnames `beta_std' = `coef_names'
    ereturn matrix beta = `beta_std'

    // 7. Save additional post-estimation results
    ereturn local cmd "nnsreg"
    ereturn local cmdline `"`nns_cmdline'"'
    ereturn local title "Nonlinear Nonparametric Statistics (NNS) regression"
    ereturn local predict "nnsreg_p"
    ereturn local depvar "`depvar'"
    ereturn local indepvar "`indepvar'"
    ereturn local method "`method'"
    ereturn local noise "`noise'"
    ereturn local partition "`partition'"
    ereturn scalar order = `order'
    ereturn scalar minobs = `minobs'
    ereturn local knots "`knots'"
    // regress supplies postestimation hooks for its original spline terms.
    // Those terms are temporary and no longer exist after nnsreg returns, so
    // advertising regress estat or margins support produces invalid calls.
    ereturn local estat_cmd "nnsreg_estat"
    ereturn local marginsok ""
    ereturn local marginsnotok "_ALL"

    // Range of the regressor over the estimation sample, so predict can flag
    // out-of-range extrapolation.
    quietly summarize `indepvar' if `touse', meanonly
    ereturn scalar xmin = r(min)
    ereturn scalar xmax = r(max)

    // Full-precision interior knot locations as a column vector.
    if `num_knots' > 0 {
        matrix colnames `knotmat' = knot
        ereturn matrix knotmat = `knotmat'
    }

    // 8. Per-segment slopes: cumulative sums of the spline coefficients, with
    // delta-method standard errors from e(V).  Same conditional-on-knots
    // caveat as the coefficient standard errors.
    tempname segmat
    mata: nns_build_segments("`touse'", "`indepvar'", "knots", "`segmat'")
    matrix colnames `segmat' = x_lo x_hi slope se n
    ereturn matrix segments = `segmat'

    // 9. Generate fitted values if requested
    if "`generate'" != "" {
        quietly predict double `generate' if `touse'
        label variable `generate' "NNS fitted values for `depvar' (`method')"
    }

    // 10. Display regression output
    local ylab : variable label `depvar'
    if `"`ylab'"' == "" local ylab "`depvar'"
    local xlab : variable label `indepvar'
    if `"`xlab'"' == "" local xlab "`indepvar'"

    di ""
    di as text "Nonlinear Nonparametric Statistics (NNS) Regression"
    di as text "Dependent variable:   " as result "`depvar'"
    di as text "Independent variable: " as result "`indepvar'"
    di as text "Method:               " as result "`method'"
    di as text "Noise reduction:      " as result "`noise'"
    di as text "Partition:            " as result "`partition'"
    di as text "Order:                " as result "`order'"
    di as text "Number of obs:        " as result "`N'"
    di as text "R-squared:            " as result %6.4f e(r2)
    di as text "Root MSE:             " as result %6.4f e(rmse)
    di ""
    ereturn display, level(`level')
    if "`method'" == "connect" {
        di as text "note: connect-fit standard errors are heuristic references from a " ///
            "rescaled OLS spline covariance; they are not estimator-specific tests."
    }

    // Optional per-segment slope table.
    if "`slopes'" != "" {
        nnsreg_slopes, level(`level')
    }

    // 11. Draw plot
    if "`noplot'" == "" {
        tempvar fitted
        quietly predict double `fitted' if `touse'

        // Texas 2036 Brand Colors:
        // Navy: "27 45 85" (RGB for #1B2D55)
        // Orange: "212 69 0" (RGB for #D44500)
        // Muted Gray: "108 122 141" (RGB for #6C7A8D)
        // Light Background: "245 247 250" (RGB for #F5F7FA)

        // Optional vertical lines at the interior knots.
        local knotopt ""
        if "`knotlines'" != "" & `num_knots' > 0 {
            local knotopt xline(`knots', lcolor("108 122 141") lpattern(shortdash) lwidth(vthin))
        }

        twoway ///
            (scatter `depvar' `indepvar' if `touse', ///
                mcolor("108 122 141") msymbol(Oh) msize(small)) ///
            (line `fitted' `indepvar' if `touse', sort ///
                lcolor("212 69 0") lwidth(medthick)), ///
            `knotopt' ///
            title("NNS regression fit, order `order'", size(medsmall)) ///
            subtitle("method(`method'), partition(`partition'), noise(`noise')", size(small)) ///
            xtitle("`xlab'", size(small)) ytitle("`ylab'", size(small)) ///
            xlabel(, labsize(small)) ylabel(, labsize(small) angle(horizontal)) ///
            legend(label(1 "Observed") label(2 "NNS fit") ring(0) position(11) ///
                rows(2) size(small) region(lcolor(none))) ///
            graphregion(color("245 247 250")) `graph_opts'
    }
end

// Display the stored e(segments) matrix as a per-segment slope table.
program define nnsreg_slopes
    syntax [, Level(cilevel)]
    if "`e(cmd)'" != "nnsreg" {
        di as error "nnsreg_slopes works only after nnsreg"
        exit 301
    }
    tempname S
    matrix `S' = e(segments)
    local nseg = rowsof(`S')
    local df_r = e(df_r)
    local tcrit = invttail(`df_r', (100 - `level') / 200)
    local nempty = 0
    local nthin = 0
    local minobs : display %9.0g e(minobs)
    local minobs = strtrim("`minobs'")

    di ""
    if "`e(method)'" == "connect" {
        di as text "Per-segment slopes (cumulative); standard errors are heuristic references"
    }
    else {
        di as text "Per-segment slopes (cumulative), conditional on the NNS knots"
    }
    di as text "{hline 73}"
    di as text %8s "segment" "  " %22s "x range" "  " %10s "slope" "  " ///
        %10s "std. err." "  " %6s "n"
    di as text "{hline 73}"
    forvalues s = 1/`nseg' {
        local lo = `S'[`s', 1]
        local hi = `S'[`s', 2]
        local sl = `S'[`s', 3]
        local se = `S'[`s', 4]
        local n  = `S'[`s', 5]
        if `n' == 0 local ++nempty
        if `n' < e(minobs) local ++nthin
        local rng : di %8.3g `lo' " to " %8.3g `hi'
        di as result %8.0f `s' "  " as text %22s "`rng'" "  " ///
            as result %10.4f `sl' "  " %10.4f `se' "  " %6.0f `n'
    }
    di as text "{hline 73}"
    if "`e(method)'" == "connect" {
        di as text "Standard errors use a rescaled OLS spline covariance and are not tests."
    }
    else {
        di as text "Standard errors are conditional on the data-derived knot locations."
    }
    if `nthin' > 0 {
        di as text "note: `nthin' reported x range(s) contain fewer than minobs(`minobs') " ///
            "observations; minobs() governs recursive " ///
            "partition groups, not the counts between adjacent knots."
    }
    if `nempty' > 0 {
        di as text "note: `nempty' reported x range(s) contain no estimation " ///
            "observations; their slopes describe interpolation between knots."
    }
end

// =============================================================================
// Mata Implementation of the NNS Core
// =============================================================================
version 16.0
mata:

// Median calculator
real scalar nns_median(real colvector v)
{
    real colvector sorted
    real scalar n, mid
    sorted = sort(v, 1)
    n = rows(sorted)
    if (n == 0) return(.)
    if (n - 2 * floor(n / 2) == 1) {
        return(sorted[(n+1)/2])
    }
    else {
        mid = n / 2
        return((sorted[mid] + sorted[mid+1]) / 2)
    }
}

// Simple Mean calculator
real scalar nns_mean(real colvector v)
{
    if (rows(v) == 0) return(.)
    return(sum(v) / rows(v))
}

// Simple Linear Regression (returns intercept, slope)
real rowvector nns_fast_lm(real colvector x, real colvector y)
{
    real scalar b, a, mx, my, var_x
    real colvector diff_x, diff_y
    real rowvector res

    mx = nns_mean(x)
    my = nns_mean(y)
    diff_x = x :- mx
    diff_y = y :- my
    var_x = sum(diff_x :* diff_x)
    if (var_x == 0) {
        b = 0
    }
    else {
        b = sum(diff_x :* diff_y) / var_x
    }
    a = my - b * mx

    res = J(1, 2, .)
    res[1] = a
    res[2] = b
    return(res)
}

// Consolidate duplicate X values by averaging Y values
void nns_consolidate(real colvector rx, real colvector ry, real colvector urx, real colvector ury)
{
    real colvector idx
    real scalar i, n, count, sum_y, p

    idx = order(rx, 1)
    rx = rx[idx]
    ry = ry[idx]
    n = rows(rx)

    if (n <= 1) {
        urx = rx
        ury = ry
        return
    }

    urx = J(n, 1, .)
    ury = J(n, 1, .)
    count = 1
    urx[1] = rx[1]
    sum_y = ry[1]
    p = 1
    for (i = 2; i <= n; i++) {
        if (rx[i] == rx[i-1]) {
            sum_y = sum_y + ry[i]
            p = p + 1
        }
        else {
            ury[count] = sum_y / p
            count = count + 1
            urx[count] = rx[i]
            sum_y = ry[i]
            p = 1
        }
    }
    ury[count] = sum_y / p

    urx = urx[1..count]
    ury = ury[1..count]
}

// Endpoint Y estimation
real scalar nns_endpoint_y(real colvector x, real colvector y, real colvector rp_x, real scalar low)
{
    real scalar boundary, reg_range, mid_range, y_boundary_val
    real colvector boundary_mask, mid_mask, y_boundary, x_mid, y_mid, boundary_y_val
    real rowvector fit_boundary, fit_mid

    if (low) {
        boundary = min(x)
        reg_range = min(rp_x)
    }
    else {
        boundary = max(x)
        reg_range = max(rp_x)
    }
    mid_range = nns_mean((boundary \ reg_range))

    if (low) {
        boundary_mask = (x :<= reg_range)
        mid_mask = (x :<= mid_range)
    }
    else {
        boundary_mask = (x :>= reg_range)
        mid_mask = (x :>= mid_range)
    }

    y_boundary = select(y, boundary_mask)
    y_mid = select(y, mid_mask)
    x_mid = select(x, mid_mask)

    if (rows(uniqrows(x_mid)) > 1 && rows(y_boundary) > 5) {
        fit_boundary = nns_fast_lm(select(x, boundary_mask), y_boundary)
        fit_mid = nns_fast_lm(x_mid, y_mid)

        y_boundary_val = fit_boundary[1] + fit_boundary[2] * boundary
        y_boundary_val = (y_boundary_val * rows(y_boundary) + (fit_mid[1] + fit_mid[2] * boundary) * rows(y_mid)) / (rows(y_boundary) + rows(y_mid))
        return(y_boundary_val)
    }

    boundary_y_val = select(y, x :== boundary)
    return(nns_mean(uniqrows(boundary_y_val)))
}

// Central point logic (returns central_x, central_y)
real rowvector nns_central_point(real colvector rp_x, real colvector rp_y, real colvector x, real colvector y)
{
    real scalar n_points, r1, r2, central_x, central_y
    real colvector row_positions, mask, x_sub, y_sub
    real rowvector res

    n_points = rows(rp_x)
    row_positions = 1::n_points
    r1 = floor(nns_median(row_positions))
    r2 = ceil(nns_median(row_positions))

    if (r1 != r2) {
        mask = (x :>= rp_x[r1]) :& (x :<= rp_x[r2])
        x_sub = select(x, mask)
        y_sub = select(y, mask)
        central_x = nns_mean((rp_x[r1] \ rp_x[r2]))
        central_y = nns_mean(y_sub)
        // No observations may lie between the two central regression points
        // when the regressor support is clustered or gapped; fall back to
        // the midpoint of the two regression points
        if (central_y >= .) central_y = nns_mean((rp_y[r1] \ rp_y[r2]))
    }
    else {
        central_x = rp_x[r1]
        central_y = rp_y[r1]
    }

    res = J(1, 2, .)
    res[1] = central_x
    res[2] = central_y
    return(res)
}

// Main NNS Regression Mata wrapper
void nns_reg_mata(string scalar depvar, string scalar indepvar, string scalar touse,
                  real scalar order, real scalar minobs, string scalar noise, string scalar partition,
                  string scalar b_name, string scalar knots_name, string scalar knotmat_name,
                  string scalar fail_name)
{
    real colvector x, y, quad_ids, unique_quads, counts, mask, x_sub, y_sub
    real scalar N, depth, i, j, K, P, cx, cy, any_split, obs_req
    real colvector rp_x, rp_y, unique_groups, rp_order
    real scalar min_x, max_x, min_y, max_y
    real rowvector central
    real colvector combined_x, combined_y, sorted_x, sorted_y
    real colvector knots, okpts
    string scalar knots_str
    real colvector rise, run, slope, beta
    real scalar alpha, M
    real matrix b_nns

    // Load data from Stata
    x = st_data(., indepvar, touse)
    y = st_data(., depvar, touse)
    N = rows(x)

    // Degenerate support: the partition needs at least a few distinct
    // regressor values to form quadrants.  Report a friendly message.
    if (rows(uniqrows(x)) < 3) {
        st_local(fail_name, "the regressor has fewer than 3 distinct values in the estimation sample")
        return
    }

    // Order 1 uses one full-sample regression point.  Each higher order adds
    // one generation of recursive partitions.
    quad_ids = J(N, 1, 1)
    obs_req = minobs // minimum group size to keep partitioning (default 8)

    // Partitioning Loop
    for (depth = 2; depth <= order; depth++) {
        unique_quads = uniqrows(quad_ids)

        counts = J(rows(unique_quads), 1, 0)
        for (i = 1; i <= rows(unique_quads); i++) {
            counts[i] = sum(quad_ids :== unique_quads[i])
        }

        any_split = 0
        for (i = 1; i <= rows(unique_quads); i++) {
            if (counts[i] > obs_req) {
                P = unique_quads[i]
                mask = (quad_ids :== P)
                x_sub = select(x, mask)

                // Compute center
                if (noise == "median") {
                    cx = nns_median(x_sub)
                    cy = nns_median(select(y, mask))
                }
                else {
                    cx = nns_mean(x_sub)
                    cy = nns_mean(select(y, mask))
                }

                // Split observations. The default follows the NNS quadrant
                // construction; xonly keeps the earlier vertical-band variant.
                for (j = 1; j <= N; j++) {
                    if (mask[j]) {
                        if (partition == "xonly") {
                            if (x[j] > cx) {
                                quad_ids[j] = 2 * P
                            }
                            else {
                                quad_ids[j] = 2 * P - 1
                            }
                        }
                        else {
                            if (x[j] <= cx && y[j] <= cy) {
                                quad_ids[j] = 4 * P - 3
                            }
                            else if (x[j] > cx && y[j] <= cy) {
                                quad_ids[j] = 4 * P - 2
                            }
                            else if (x[j] <= cx && y[j] > cy) {
                                quad_ids[j] = 4 * P - 1
                            }
                            else {
                                quad_ids[j] = 4 * P
                            }
                        }
                    }
                }
                any_split = 1
            }
        }
        if (!any_split) break
    }

    // Group by the final partition IDs to compute regression points
    unique_groups = uniqrows(quad_ids)
    K = rows(unique_groups)
    rp_x = J(K, 1, .)
    rp_y = J(K, 1, .)
    for (i = 1; i <= K; i++) {
        mask = (quad_ids :== unique_groups[i])
        x_sub = select(x, mask)
        y_sub = select(y, mask)
        if (noise == "median") {
            rp_x[i] = nns_median(x_sub)
            rp_y[i] = nns_median(y_sub)
        }
        else {
            rp_x[i] = nns_mean(x_sub)
            rp_y[i] = nns_mean(y_sub)
        }
    }

    rp_order = order(rp_x, 1)
    rp_x = rp_x[rp_order]
    rp_y = rp_y[rp_order]

    // Calculate endpoints
    min_x = min(x)
    max_x = max(x)
    min_y = nns_endpoint_y(x, y, rp_x, 1)
    max_y = nns_endpoint_y(x, y, rp_x, 0)

    // Calculate central point
    central = nns_central_point(rp_x, rp_y, x, y)

    // Combine all regression points
    combined_x = (min_x \ max_x \ central[1] \ rp_x)
    combined_y = (min_y \ max_y \ central[2] \ rp_y)

    // Guard: drop any candidate regression point with a missing coordinate
    // so a degenerate partition cannot propagate missing values into e(b)
    okpts = (combined_x :< .) :& (combined_y :< .)
    combined_x = select(combined_x, okpts)
    combined_y = select(combined_y, okpts)

    // Consolidate duplicate X values
    sorted_x = J(0, 1, .)
    sorted_y = J(0, 1, .)
    nns_consolidate(combined_x, combined_y, sorted_x, sorted_y)

    // If consolidation collapses to a single fit point, the segment
    // construction below has no run to work with; report cleanly.
    if (rows(sorted_x) < 2) {
        st_local(fail_name, "the partition produced fewer than two distinct fit points; the regressor may be nearly constant")
        return
    }

    // The interior knots are sorted_x[2..M+1] where M = rows(sorted_x) - 2
    M = rows(sorted_x) - 2
    if (M > 0) {
        knots = sorted_x[2..M+1]
        knots_str = ""
        for (i = 1; i <= M; i++) {
            knots_str = knots_str + " " + sprintf("%20.15f", knots[i])
        }
        st_local(knots_name, strtrim(knots_str))
        st_matrix(knotmat_name, knots)
    }
    else {
        st_local(knots_name, "")
    }

    // Compute coefficients for connect method
    // slopes are dy/dx
    rise = sorted_y[2..rows(sorted_y)] - sorted_y[1..rows(sorted_y)-1]
    run = sorted_x[2..rows(sorted_x)] - sorted_x[1..rows(sorted_x)-1]
    slope = J(rows(run), 1, 0)
    for (i = 1; i <= rows(run); i++) {
        if (run[i] == 0) {
            slope[i] = 0
        }
        else {
            slope[i] = rise[i] / run[i]
        }
    }

    // constant alpha
    alpha = sorted_y[1] - slope[1] * sorted_x[1]

    // beta
    beta = J(rows(slope), 1, 0)
    beta[1] = slope[1]
    for (i = 2; i <= rows(slope); i++) {
        beta[i] = slope[i] - slope[i-1]
    }

    // construct b_nns matrix (row vector)
    b_nns = (beta' , alpha)
    st_matrix(b_name, b_nns)
}

// Calculate NNS fitted values (Mata helper)
void nns_calc_fitted(string scalar sp_vars, string scalar b_name, string scalar touse, string scalar y_hat_var)
{
    real matrix X_spline
    real rowvector b
    real colvector y_hat
    string rowvector vars

    vars = tokens(sp_vars)
    X_spline = st_data(., vars, touse)
    X_spline = X_spline , J(rows(X_spline), 1, 1) // add constant column

    b = st_matrix(b_name)
    y_hat = X_spline * b'

    st_store(., y_hat_var, touse, y_hat)
}

// Build the per-segment slope matrix: [x_lo, x_hi, slope, se, n], one row per
// segment.  Slope in segment s is the cumulative sum of the first s spline
// coefficients (first-segment slope plus the changes at earlier knots); its
// standard error is a' V a for the cumulative selector a, the matrix analogue
// of the lincom the article performs by hand.
void nns_build_segments(string scalar touse, string scalar xvar, string scalar knotloc, string scalar outname)
{
    real colvector x, knots
    real matrix b, V, out
    real scalar M, s, i, xmin, xmax, lo, hi
    real rowvector a
    string scalar kl

    x = st_data(., xvar, touse)
    xmin = min(x)
    xmax = max(x)
    b = st_matrix("e(b)")   // 1 x (M+2): [indepvar, knot1..knotM, _cons]
    V = st_matrix("e(V)")   // (M+2) x (M+2)

    kl = st_local(knotloc)
    if (strtrim(kl) == "") knots = J(0, 1, .)
    else knots = strtoreal(tokens(kl))'
    M = rows(knots)

    out = J(M + 1, 5, .)
    for (s = 1; s <= M + 1; s++) {
        if (s == 1) lo = xmin
        else lo = knots[s - 1]
        if (s == M + 1) hi = xmax
        else hi = knots[s]

        a = J(1, cols(b), 0)
        for (i = 1; i <= s; i++) a[i] = 1

        out[s, 1] = lo
        out[s, 2] = hi
        out[s, 3] = (a * b')[1, 1]
        out[s, 4] = sqrt((a * V * a')[1, 1])
        if (s < M + 1) out[s, 5] = sum((x :>= lo) :& (x :< hi))
        else out[s, 5] = sum((x :>= lo) :& (x :<= hi))
    }
    st_matrix(outname, out)
}

end
