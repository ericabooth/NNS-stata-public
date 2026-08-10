clear all
macro drop _all
set more off
set seed 20360702
adopath + "code"

tempname bench
global BENCH `bench'
postfile `bench' str24 scenario str32 model double r2 rmse_y rmse_true using "work/benchmark_results.dta", replace

capture program drop _metrics
program define _metrics
    syntax, Scenario(string) Model(string) Fitvar(name) Truevar(name) Yvar(name)
    tempvar resy rest tssy
    quietly gen double `resy' = (`yvar' - `fitvar')^2
    quietly gen double `rest' = (`truevar' - `fitvar')^2
    quietly summarize `resy', meanonly
    local rmse_y = sqrt(r(mean))
    quietly summarize `rest', meanonly
    local rmse_true = sqrt(r(mean))
    quietly summarize `yvar', meanonly
    local ybar = r(mean)
    quietly gen double `tssy' = (`yvar' - `ybar')^2
    quietly summarize `resy', meanonly
    local rss = r(sum)
    quietly summarize `tssy', meanonly
    local tss = r(sum)
    local r2 = 1 - (`rss' / `tss')
    post $BENCH ("`scenario'") ("`model'") (`r2') (`rmse_y') (`rmse_true')
end

capture program drop _run_models
program define _run_models
    syntax, Scenario(string) Xvar(name) Yvar(name) Truevar(name)
    regress `yvar' `xvar'
    predict double fit_linear, xb
    _metrics, scenario("`scenario'") model("Linear regression") fitvar(fit_linear) truevar(`truevar') yvar(`yvar')

    quietly npregress kernel `yvar' `xvar'
    predict double fit_np, mean
    _metrics, scenario("`scenario'") model("npregress kernel") fitvar(fit_np) truevar(`truevar') yvar(`yvar')

    foreach p in xy xonly {
        foreach o in 2 3 4 5 {
            capture noisily nnsreg `yvar' `xvar', order(`o') method(connect) partition(`p') noplot generate(fit_conn_`p'_`o')
            if !_rc {
                _metrics, scenario("`scenario'") model("NNS connect `p' o`o'") fitvar(fit_conn_`p'_`o') truevar(`truevar') yvar(`yvar')
            }
            capture noisily nnsreg `yvar' `xvar', order(`o') method(ols) partition(`p') noplot generate(fit_ols_`p'_`o')
            if !_rc {
                _metrics, scenario("`scenario'") model("NNS OLS `p' o`o'") fitvar(fit_ols_`p'_`o') truevar(`truevar') yvar(`yvar')
            }
        }
    }
end

// Smooth S-curve with moderate noise
clear
set obs 250
gen x = 5 + 15 * runiform()
gen true = 40 + 45 / (1 + exp(-0.8 * (x - 11)))
gen y = true + rnormal(0, 4)
_run_models, scenario("s_curve_noise") xvar(x) yvar(y) truevar(true)

// Threshold with moderate noise
clear
set obs 300
gen x = 80 * runiform()
gen true = cond(x <= 35, 8, 8 + 0.4 * (x - 35) + 0.01 * (x - 35)^2)
gen y = true + rnormal(0, 2.5)
_run_models, scenario("threshold_noise") xvar(x) yvar(y) truevar(true)

// High signal sine on a deterministic grid
clear
set obs 300
gen x = 4 * _pi * (_n - 1) / (_N - 1)
gen true = sin(x)
gen y = true + rnormal(0, .04)
_run_models, scenario("sine_grid_low_noise") xvar(x) yvar(y) truevar(true)

// High signal stochastic sine
clear
set obs 300
gen x = 4 * _pi * runiform()
gen true = sin(x)
gen y = true + rnormal(0, .04)
_run_models, scenario("sine_random_low_noise") xvar(x) yvar(y) truevar(true)

// Clustered sharp dip inspired by NNS segmentation examples
clear
set obs 260
gen mix = runiform()
gen x = cond(mix < .45, rnormal(22, 3), cond(mix < .75, rnormal(36, 2.2), rnormal(58, 4)))
replace x = max(5, min(70, x))
gen true = 18 + .18*x + 10*exp(-((x-24)^2)/35) - 16*exp(-((x-35)^2)/18) + 8*exp(-((x-58)^2)/50)
gen y = true + rnormal(0, 1.2)
_run_models, scenario("clustered_peak_dip") xvar(x) yvar(y) truevar(true)

postclose `bench'
use "work/benchmark_results.dta", clear
format r2 rmse_y rmse_true %8.4f
sort scenario rmse_true
by scenario: gen rank_true = _n
by scenario: gen rank_y = _n
list if rank_true <= 6, sepby(scenario) abbreviate(28)
export delimited using "work/benchmark_results.csv", replace
