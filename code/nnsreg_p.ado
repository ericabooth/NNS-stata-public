*! version 1.1.0  17jul2026
*! nnsreg_p: Post-estimation predict program for nnsreg
*! Author: Eric A. Booth (Senior Researcher, Texas 2036, eric.a.booth@gmail.com)
*! 1.1.0: note out-of-range scoring, since the terminal spline segments
*!        extrapolate linearly beyond the estimation-sample range of the
*!        regressor.

program define nnsreg_p
    version 16.0
    
    // Check if called after nnsreg
    if "`e(cmd)'" != "nnsreg" {
        error 301
    }
    
    syntax newvarname [if] [in] [, xb residuals]
    
    // Default is xb (fitted values)
    if "`xb'" == "" & "`residuals'" == "" {
        local xb "xb"
    }
    if "`xb'" != "" & "`residuals'" != "" {
        di as error "options xb and residuals are mutually exclusive"
        exit 198
    }

    // Retrieve model info from e()
    local knots "`e(knots)'"
    local indepvar "`e(indepvar)'"
    local depvar "`e(depvar)'"
    
    // Mark sample
    marksample touse, novarlist
    quietly replace `touse' = 0 if missing(`indepvar')

    // Recreate the spline basis locally
    local num_knots : word count `knots'
    
    tempvar sp_1
    quietly gen double `sp_1' = `indepvar' if `touse'
    local sp_vars `sp_1'
    
    forvalues j = 1/`num_knots' {
        local k : word `j' of `knots'
        tempvar sp_`=`j'+1'
        quietly gen double `sp_`=`j'+1'' = cond(`indepvar' > `k', `indepvar' - `k', 0) if `touse'
        local sp_vars `sp_vars' `sp_`=`j'+1''
    }

    // Score the prediction using e(b)
    tempname b
    matrix `b' = e(b)
    
    // Rename columns of b to match our new temporary spline variables
    matrix colnames `b' = `sp_vars' _cons
    
    // Note out-of-range scoring: the terminal segments extrapolate linearly
    // beyond the estimation-sample range stored in e(xmin)/e(xmax).
    if "`e(xmin)'" != "" & "`e(xmax)'" != "" {
        quietly count if `touse' & (`indepvar' < `e(xmin)' | `indepvar' > `e(xmax)')
        if r(N) > 0 {
            di as text "note: `r(N)' observation(s) scored outside the " ///
                "estimation-sample range of `indepvar'; the terminal " ///
                "segments extrapolate linearly."
        }
    }

    if "`xb'" != "" {
        // Fitted values
        matrix score `typlist' `varlist' = `b' if `touse'
        label variable `varlist' "NNS fitted values for `depvar'"
    }
    else if "`residuals'" != "" {
        // Residuals
        tempvar yhat
        matrix score double `yhat' = `b' if `touse'
        gen `typlist' `varlist' = `depvar' - `yhat' if `touse'
        label variable `varlist' "NNS residuals for `depvar'"
    }
end
