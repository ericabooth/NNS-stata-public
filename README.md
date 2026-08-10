# nnsreg

`nnsreg` is a Stata implementation of univariate Nonlinear Nonparametric Statistics (NNS) regression. It partitions the data, computes regression points, and represents the fitted curve as a linear spline with data-driven knots.

The package is designed for applied users who want to inspect where a fitted relationship changes, identify relatively flat or steep ranges, and check whether those features survive sensitivity checks.

## What Questions Does It Help Answer?

- Where does the fitted relationship change direction or slope?
- Which parts of the observed support look flat, steep, or threshold-like?
- Are the slope changes robust to `order()`, `noise(mean|median)`, and `partition(xy|xonly)`?
- Do segmented fits recover a known or expected pattern better than a linear model or kernel smoother?
- Can the shape be reported with stored estimates, predicted values, residuals, and coefficient plots?

This makes `nnsreg` useful when the analyst needs more than a smooth curve: for example, clustered support, sharp local features, threshold-style patterns, or a reportable set of slope-change components.

## Install

From Stata:

```stata
net install nnsreg, from("https://raw.githubusercontent.com/ericbooth/NNS-stata-public/main/github") replace
```

For local development from this folder:

```stata
adopath + "code"
which nnsreg
help nnsreg
```

## Syntax

```stata
nnsreg depvar indepvar [if] [in] [, order(integer) noise(mean|median) ///
    method(connect|ols) partition(xy|xonly) minobs(integer) level(cilevel) ///
    vce(robust|cluster clustvar) slopes noplot knotlines graph_opts(string) ///
    generate(newvar)]
```

Key options:

- `order(#)` sets the partition depth. Larger values allow more segments.
- `noise(mean)` uses means to form partition centers and regression points.
- `noise(median)` is a robustness check when outliers or skewness may matter.
- `method(connect)` connects NNS regression points directly.
- `method(ols)` estimates an OLS linear spline using NNS-derived knots.
- `partition(xy)` is the default NNS-style quadrant partition.
- `partition(xonly)` forms vertical bands and is useful as a sensitivity check.
- `minobs(#)` sets the smallest group the partition keeps splitting (default 8). An adjacent-knot range can still contain fewer observations, or none.
- `vce(robust|cluster clustvar)` requests robust or cluster-robust standard errors for the least-squares stage (`method(ols)` only).
- `slopes` reports per-segment slopes with delta-method standard errors (also stored in `e(segments)`).
- `knotlines` overlays the interior knots on the default plot.
- `generate(newvar)` stores fitted values for the estimation sample.

## Example Workflow

Run the full simulation, comparison, and figure script:

```stata
do examples/test_nnsreg.do
```

The script:

- Simulates an education example with diminishing returns.
- Simulates a water infrastructure example with threshold behavior.
- Simulates a clustered peak/dip benchmark inspired by NNS segmentation examples.
- Compares `regress`, `npregress kernel`, `nnsreg, method(connect)`, and `nnsreg, method(ols)`.
- Runs NNS sensitivity checks for `order()`, `noise()`, and partition choices.
- Exports single-method fit panels, a coefficient plot with knot-located labels, and a knot-location figure to `manuscript/`.
- Repeats the whole comparison across 100 Monte Carlo replications of each design and writes `manuscript/table_mc_summary.csv`.
- Writes single-run diagnostics to `manuscript/table_model_diagnostics.csv`.

## Interpreting The Examples

The examples are simulated. They test whether the commands recover known nonlinear shapes and do not estimate real Texas policy effects.

The 100-replication Monte Carlo study (seeds recorded in the script) compares `nnsreg` with `npregress kernel` and `npregress series`:

- In the smooth education and water designs, both `npregress` estimators recover the true curve more accurately on average; NNS OLS is within about 0.06 to 0.14 RMSE.
- In the clustered peak/dip design, `npregress series` reaches the lowest average error but has the widest distribution (standard deviation 0.33 across draws, larger than the kernel's); NNS OLS nearly matches it on average with about one-fifth to one-quarter of the standard deviation.
- Higher `order()` helps only where the true curve has sharp local features; in the smooth designs it tracks noise.

On the motorcycle data, repeated five-fold cross-validation gives the lowest mean held-out RMSE to `npregress series` (24.08), with `npregress kernel` (24.52) and NNS OLS (24.55) close behind. The NNS output additionally reports range slopes and support counts for descriptive follow-up.

That is the intended use case: consider `nnsreg` when segmented structure or slope-change interpretation matters, particularly as a differently tuned comparison on clustered support. It reports a compact segment summary with named ranges, slopes, and support counts. Use validation and comparison fits to decide whether it adds value in a given application.

## Postestimation

After fitting:

```stata
nnsreg staar_pass spend_pupil, order(3) method(ols) partition(xy) generate(yhat)
estimates store nns_ols
predict resid, residuals
estimates table nns_ols, b se stats(r2 rmse N)
```

The first coefficient is the first segment slope. Knot coefficients are changes in slope after adaptive knots.

## Repository Layout

- `code/`: Stata package files (`nnsreg.ado`, `nnsreg_p.ado`,
  `nnsreg_estat.ado`, `nnsreg.sthlp`).
- `examples/`: simulated data generator, full verification script, the
  `motorcycle_example.do` real-data illustration (repeated k-fold
  cross-validation), `robustness.do` (interval coverage and an alternative
  clustered design), and `r_fidelity.R` (comparison against the R NNS package).
- `data/`: simulated example datasets.
- `manuscript/`: Stata Journal draft, figures, and submission packet.
- `webpage/`: webdoc2 page and statashiny companion.
- `github/`: `net install` files.
- `literature/`: background literature organized into `core_nns/` and `comparators/`.

## Caveats

The current package is univariate. The `method(ols)` standard errors condition on knots selected from the same data. The `method(connect)` intervals are heuristic references based on a rescaled OLS covariance matrix, not estimator-specific tests. A bend in the fitted line is descriptive unless paired with a research design that supports a causal claim. Empty or thinly supported ranges should be treated as prompts for inspection, not as local evidence.

When prediction is the main goal, compare `nnsreg` against `npregress kernel`, splines, and simpler baselines using a common validation metric. A future machine-learning version would be most useful if it adds validation-based tuning, multivariate partitioning, and uncertainty summaries.

## Citation
Booth, Eric A. “Reportable Segment Slopes for Nonlinear Relationships in Stata (nnsreg)”. MetaArXiv, pre-print (Under Review at Stata Journal)  osf.io/preprints/metaarxiv/dftr7_v1.

## Author

Eric A. Booth  
Senior Researcher, Texas 2036  
eric.a.booth@gmail.com
