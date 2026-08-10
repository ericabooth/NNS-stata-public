nnsreg: Reportable segment slopes for nonlinear relationships in Stata
Stata Journal submission packet
Author: Eric A. Booth (eric.a.booth@gmail.com)

This packet is flat: every file sits in this directory, nnsreg.tex compiles in
place, and every command in the manuscript reproduces as typed from here.


CONTENTS
--------

Manuscript and Stata Journal LaTeX files
  nnsreg.tex             manuscript source (single self-contained file)
  nnsreg.bib             bibliography
  nnsreg.pdf             compiled manuscript
  statapress.cls, sj.sty, stata.sty, pagedims.sty, sj.bst
                         Stata Journal document class and style files

Package files (the command)
  nnsreg.ado             the command
  nnsreg_p.ado           the predict program
  nnsreg_estat.ado       dispatcher that rejects unsupported estat calls
  nnsreg.sthlp           help file

Reproducibility do-files (path-flattened for this directory)
  simulated_data.do      defines nns_simdata and writes the three example
                         datasets with fixed seeds
  test_nnsreg.do         verification tests, the three worked examples, the
                         100-replication Monte Carlo study, all simulated-data
                         figures, and the diagnostic tables
  make_sjlogs.do         regenerates the four .log.tex output excerpts
  motorcycle_example.do  the real-data illustration: fit, repeated
                         cross-validation, the segment-support figure, and
                         the motorcycle tables
  robustness.do          empirical coverage of the conditional intervals and
                         an alternative higher-noise clustered design

Example datasets (regenerable from simulated_data.do)
  texas_education.dta, texas_water.dta, clustered_peak_dip.dta

Stata output excerpts included by the manuscript
  nnsreg_ex_edu.log.tex, nnsreg_ex_plateau.log.tex,
  nnsreg_ex_water.log.tex, nnsreg_ex_slopes.log.tex

Figures included by the manuscript (PDF, grayscale)
  fig_edu_fits.pdf, fig_edu_coef.pdf, fig_edu_sensitivity.pdf,
  fig_water_fits.pdf, fig_water_knots.pdf, fig_water_sensitivity.pdf,
  fig_cluster_fits.pdf, fig_cluster_sensitivity.pdf, fig_mc_rmse.pdf,
  fig_moto_fits.pdf, fig_moto_knots.pdf

Editorial correspondence
  cover_letter.txt       letter to the editors


BUILD THE PDF
-------------

  latexmk -pdf nnsreg.tex


REPRODUCE THE RESULTS
---------------------

From Stata, in this directory:

  do simulated_data.do        // write the three example datasets
  do test_nnsreg.do           // tests, examples, Monte Carlo, figures, tables
  do make_sjlogs.do           // the four .log.tex output excerpts
  do motorcycle_example.do    // real-data fit, CV table, segment figure
  do robustness.do            // interval coverage + alternative clustered design

test_nnsreg.do also calls simulated_data.do, installs coefplot from SSC if it
is not already present, and uses the built-in npregress kernel and npregress
series. The Monte Carlo study runs 100 replications per design (about three
minutes). motorcycle_example.do loads the motorcycle data with webuse and runs
a repeated five-fold cross-validation. Every figure and table regenerates from
these scripts with the seeds set inside them. nnsreg requires Stata 16 or
later.


NOTES
-----

The simulated examples use known generating curves so that each method can be
scored against the true curve. The variable names borrow from Texas policy
settings for concreteness; the values are synthetic and stylized, and no claim
is made about actual districts, assessments, or water systems. The real-data
illustration uses the motorcycle-impact benchmark distributed by Stata Press
(webuse motorcycle).
