clear all
set more off
adopath ++ "code"
set seed 1729
set obs 250
generate double x = runiform()
generate double y = sin(6*x) + rnormal()
generate long cid = ceil(_n/10)
replace cid = . in 241/250

capture noisily nnsreg y x, method(ols) vce(hc2) noplot
assert _rc == 198
capture noisily nnsreg y x, method(ols) vce(bootstrap) noplot
assert _rc == 198
capture noisily nnsreg y x, method(ols) vce(jackknife) noplot
assert _rc == 198
capture noisily nnsreg y x, method(ols) vce(robust) noplot
assert _rc == 0
capture noisily nnsreg y x, method(ols) vce(cluster cid) generate(yhat) noplot
assert _rc == 0
assert e(N) == 240
count if e(sample)
assert r(N) == 240
count if !missing(yhat)
assert r(N) == 240
matrix S = e(segments)
mata: st_numscalar("seg_n", sum(st_matrix("S")[,5]))
assert scalar(seg_n) == 240
capture noisily nnsreg y x, method(ols) slopes noplot
assert _rc == 0
display "VCE FIX PROBE PASSED"
