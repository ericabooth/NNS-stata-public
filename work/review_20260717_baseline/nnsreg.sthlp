{smcl}
{* *! version 1.3.0  17jul2026}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[R] npregress kernel" "help npregress kernel"}{...}
{vieweralsosee "[R] lpoly" "help lpoly"}{...}
{vieweralsosee "[R] mkspline" "help mkspline"}{...}
{vieweralsosee "coefplot" "help coefplot"}{...}
{viewerjumpto "Syntax" "nnsreg##syntax"}{...}
{viewerjumpto "Description" "nnsreg##description"}{...}
{viewerjumpto "Options" "nnsreg##options"}{...}
{viewerjumpto "Examples" "nnsreg##examples"}{...}
{viewerjumpto "Stored results" "nnsreg##results"}{...}
{viewerjumpto "Cautions" "nnsreg##cautions"}{...}
{viewerjumpto "Author" "nnsreg##author"}{...}

{title:Title}

{phang}
{bf:nnsreg} {hline 2} Nonlinear Nonparametric Statistics (NNS) regression in Stata

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:nnsreg}
{depvar}
{indepvar}
{ifin}
[{cmd:,} {it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt ord:er(#)}}hierarchical partition depth; default is {cmd:order(3)}{p_end}
{synopt:{opt noi:se(string)}}partition center and regression-point summary: {cmd:mean} or {cmd:median}; default is {cmd:noise(mean)}{p_end}
{synopt:{opt meth:od(string)}}fit type: {cmd:connect} or {cmd:ols}; default is {cmd:method(connect)}{p_end}
{synopt:{opt part:ition(string)}}partition style: {cmd:xy} or {cmd:xonly}; default is {cmd:partition(xy)}{p_end}
{synopt:{opt min:obs(#)}}minimum group size that a partition keeps splitting; default is {cmd:minobs(8)}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt vce(vcetype)}}{it:vcetype} may be {opt r:obust} or {opt cl:uster} {it:clustvar}; allowed only with {cmd:method(ols)}{p_end}
{synopt:{opt slop:es}}report per-segment slopes and their standard errors{p_end}

{syntab:Plot}
{synopt:{opt nop:lot}}suppress the default scatter-and-fit graph{p_end}
{synopt:{opt knotl:ines}}overlay vertical lines at the interior knots on the default graph{p_end}
{synopt:{opt graph_opts(string)}}pass additional {help twoway} graph options{p_end}

{syntab:Generate}
{synopt:{opt gen:erate(newvar)}}create a fitted-value variable{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
After estimation, {cmd:predict} {newvar} [{cmd:,} {cmd:xb} {cmd:residuals}]
computes fitted values ({cmd:xb}, the default) or residuals.  The spline
basis is rebuilt from {cmd:e(knots)}, so {cmd:predict} also scores
observations outside the estimation sample.  {cmd:nnsreg} requires Stata 16
or later.

{marker description}{...}
{title:Description}

{pstd}
{cmd:nnsreg} fits a univariate Nonlinear Nonparametric Statistics (NNS) regression curve.
The command recursively partitions the data using partial-moment quadrants, computes a regression point inside each partition, and represents the fitted curve as a linear spline with data-derived knots.

{pstd}
The command is most useful when an analyst wants to describe where a relationship bends, which ranges are relatively flat or steep, and whether those features survive simple option checks.
In simulation comparisons distributed with the package, {helpb npregress kernel} and {cmd:npregress series} recover smooth conditional means more accurately on average, while {cmd:nnsreg}'s {cmd:method(ols)} variant matches the most accurate smoother on clustered, sharply featured data and is far more stable across replications.
{cmd:nnsreg} also fits without a bandwidth choice and yields a compact, testable segment summary that {cmd:npregress} does not report.

{pstd}
The coefficients are posted as a spline basis.  The first coefficient is the first-segment slope.
Coefficients named {cmd:knot1}, {cmd:knot2}, and so on are changes in slope after the NNS-derived knots, whose locations are stored in {cmd:e(knots)}.
This design supports {helpb estimates store}, {helpb estimates table}, {helpb predict}, {helpb lincom}, and coefficient plotting workflows.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt order(#)} specifies the partition depth.  A larger value allows more line segments.  This can follow curvature more closely, but it can also make the fit more sensitive to local noise.  The default is {cmd:order(3)}.

{phang}
{opt noise(string)} specifies the center used for partitioning and regression points.  {cmd:noise(mean)}, the default, uses arithmetic means.  {cmd:noise(median)} uses medians and is a useful sensitivity check when the regressor or outcome is skewed or contains outlying values.

{phang}
{opt method(string)} specifies the fit:

{phang2}
{cmd:method(connect)} connects NNS regression points.  This is the default and is closest to the line-segment construction described in the NNS literature.

{phang2}
{cmd:method(ols)} uses the same NNS-derived knots but estimates the linear spline coefficients by ordinary least squares.  This typically improves in-sample fit because coefficients are optimized after the knots are chosen.

{phang}
{opt partition(string)} specifies the partition style.  {cmd:partition(xy)}, the default, splits eligible groups using both the regressor and the outcome, following the NNS quadrant construction.  {cmd:partition(xonly)} splits only along the regressor and can be useful as a vertical-band sensitivity check.

{phang}
{opt minobs(#)} sets the smallest group the recursive partition keeps splitting; a group is split only if it holds more than {it:#} observations.  The default is {cmd:minobs(8)}.  Larger values produce fewer, longer segments and are a useful sensitivity check alongside {opt order()}.

{dlgtab:Reporting}

{phang}
{opt level(#)} sets the confidence level for the reported intervals; the default is {cmd:level(95)} or as set by {helpb set level}.

{phang}
{opt vce(vcetype)} passes a variance estimator through to the least-squares stage.  {it:vcetype} may be {opt robust} or {opt cluster} {it:clustvar}.  It is allowed only with {cmd:method(ols)}, whose coefficients are ordinary least squares on the spline basis; the standard errors remain conditional on the data-derived knots.

{phang}
{opt slopes} reports a table of per-segment slopes, each the cumulative sum of the spline coefficients up to that segment, with delta-method standard errors.  These are the same quantities the manual {helpb lincom} recipe below produces; the table is also stored in {cmd:e(segments)}.

{dlgtab:Plot}

{phang}
{opt noplot} suppresses the default graph.

{phang}
{opt knotlines} overlays vertical dashed lines at the interior knot locations on the default graph.

{phang}
{opt graph_opts(string)} passes standard {help twoway} options to the default graph.

{dlgtab:Generate}

{phang}
{opt generate(newvar)} creates a fitted-value variable.

{marker examples}{...}
{title:Examples}

{pstd}The examples use the simulated datasets distributed with the package
(created by {cmd:simulated_data.do} with fixed seeds).{p_end}

{pstd}Education example: fit NNS connect and NNS OLS, then compare{p_end}
{phang2}{stata "use texas_education, clear" : . use texas_education, clear}{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) method(connect) generate(fit_connect)" : . nnsreg staar_pass spend_pupil, order(3) method(connect) generate(fit_connect)}{p_end}
{phang2}{stata "estimates store edu_connect" : . estimates store edu_connect}{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) method(ols) generate(fit_ols)" : . nnsreg staar_pass spend_pupil, order(3) method(ols) generate(fit_ols)}{p_end}
{phang2}{stata "estimates store edu_ols" : . estimates store edu_ols}{p_end}
{phang2}{stata "estimates table edu_connect edu_ols, b se stats(r2 rmse N)" : . estimates table edu_connect edu_ols, b se stats(r2 rmse N)}{p_end}

{pstd}Test the slope of the final segment (a plateau test): the terminal
slope is the first-segment coefficient plus all knot coefficients{p_end}
{phang2}{stata "estimates restore edu_ols" : . estimates restore edu_ols}{p_end}
{phang2}{stata "lincom spend_pupil + knot1 + knot2 + knot3 + knot4 + knot5 + knot6 + knot7 + knot8 + knot9" : . lincom spend_pupil + knot1 + knot2 + ... }{p_end}

{pstd}Report the per-segment slopes directly (same quantities as the lincom
recipe, read from {cmd:e(segments)}){p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) method(ols) slopes" : . nnsreg staar_pass spend_pupil, order(3) method(ols) slopes}{p_end}

{pstd}Cluster-robust standard errors for the least-squares fit{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) method(ols) vce(cluster district_id)" : . nnsreg staar_pass spend_pupil, order(3) method(ols) vce(cluster district_id)}{p_end}

{pstd}Sensitivity checks: change order, center, partition style, and minimum group size{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(2) noplot" : . nnsreg staar_pass spend_pupil, order(2) noplot}{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) noise(median) noplot" : . nnsreg staar_pass spend_pupil, order(3) noise(median) noplot}{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) method(ols) partition(xonly) noplot" : . nnsreg staar_pass spend_pupil, order(3) method(ols) partition(xonly) noplot}{p_end}
{phang2}{stata "nnsreg staar_pass spend_pupil, order(3) minobs(20) noplot" : . nnsreg staar_pass spend_pupil, order(3) minobs(20) noplot}{p_end}

{pstd}Compare with linear regression and kernel regression{p_end}
{phang2}{stata "regress staar_pass spend_pupil" : . regress staar_pass spend_pupil}{p_end}
{phang2}{stata "npregress kernel staar_pass spend_pupil" : . npregress kernel staar_pass spend_pupil}{p_end}

{pstd}Prediction after {cmd:nnsreg}{p_end}
{phang2}{stata "predict yhat" : . predict yhat}{p_end}
{phang2}{stata "predict resid, residuals" : . predict resid, residuals}{p_end}

{pstd}Run the full verification, comparison, and Monte Carlo script{p_end}
{phang2}{stata "do test_nnsreg.do" : . do test_nnsreg.do}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:nnsreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(r2)}}coefficient of determination, 1 - sum(e^2)/sum((y-ybar)^2){p_end}
{synopt:{cmd:e(rmse)}}root mean squared error, sqrt(sum(e^2)/df_r){p_end}
{synopt:{cmd:e(order)}}partition depth{p_end}
{synopt:{cmd:e(minobs)}}minimum group size for splitting{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(xmin)}}minimum of the regressor in the estimation sample{p_end}
{synopt:{cmd:e(xmax)}}maximum of the regressor in the estimation sample{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:nnsreg}{p_end}
{synopt:{cmd:e(predict)}}{cmd:nnsreg_p}{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(indepvar)}}independent variable name{p_end}
{synopt:{cmd:e(method)}}{cmd:connect} or {cmd:ols}{p_end}
{synopt:{cmd:e(noise)}}{cmd:mean} or {cmd:median}{p_end}
{synopt:{cmd:e(partition)}}{cmd:xy} or {cmd:xonly}{p_end}
{synopt:{cmd:e(knots)}}space-separated interior knots{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(knotmat)}}interior knot locations (full precision, column vector){p_end}
{synopt:{cmd:e(segments)}}per-segment table: x_lo, x_hi, slope, se, n, one row per segment{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
For {cmd:method(connect)}, {cmd:e(r2)} and {cmd:e(rmse)} are computed from
raw (uncentered) squared residuals because the connected fit is not least
squares, so its residuals need not average zero.

{marker cautions}{...}
{title:Cautions}

{pstd}
{cmd:nnsreg} is descriptive unless paired with a research design.  A bend or threshold in the fitted curve is not a causal estimate by itself.

{pstd}
Standard errors are conditional on the data-derived knots, and for {cmd:method(connect)} they come from a rescaled OLS spline covariance matrix.  Read the intervals as descriptive summaries of segment stability rather than exact inference.  Use graphs, comparison models, and option sensitivity checks before making substantive claims.  When prediction accuracy is central, compare {cmd:nnsreg} with {helpb npregress kernel}, splines, and simpler baselines using a common validation metric.

{pstd}
{helpb margins} is not supported: the coefficient names ({cmd:knot1}, {cmd:knot2}, and so on) are spline-basis terms, not variables in the dataset.  Use {helpb lincom} on the stored coefficients, or the {opt slopes} option and {cmd:e(segments)}, to report marginal slopes.  {cmd:predict} rebuilds the spline basis from {cmd:e(knots)} and scores observations outside the estimation sample; because the terminal segments extrapolate linearly, {cmd:predict} notes when it scores points beyond {cmd:e(xmin)} and {cmd:e(xmax)}.

{marker author}{...}
{title:Author}

{pstd}
Eric A. Booth{break}
Texas 2036{break}
eric.a.booth@gmail.com{break}
https://github.com/ericbooth/NNS-stata-public

{title:References}

{p 4 8 2}
Hayfield, T., and J. S. Racine. 2008. Nonparametric econometrics: The np package. {it:Journal of Statistical Software} 27(5): 1-32.

{p 4 8 2}
Jann, B. 2014. Plotting regression coefficients and other estimates. {it:Stata Journal} 14(4): 708-737.

{p 4 8 2}
Vinod, H. D., and F. Viole. 2018. Nonparametric regression using clusters. {it:Computational Economics} 52(4): 1317-1334.

{p 4 8 2}
Viole, F. 2026. NNS: Nonlinear Nonparametric Statistics. R package version 13.0. https://CRAN.R-project.org/package=NNS

{p 4 8 2}
Viole, F., and D. Nawrocki. 2013. {it:Nonlinear Nonparametric Statistics: Using Partial Moments}. CreateSpace Independent Publishing Platform.
{p_end}
