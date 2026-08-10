*! version 1.3.0  05jul2026
*! Simulated data generator for the nnsreg Stata Journal manuscript.
*! Defines nns_simdata, a reusable data-generating program, and writes the
*! three canonical example datasets to examples/.  The Monte Carlo study in
*! examples/test_nnsreg.do calls nns_simdata with fresh seeds, so the single
*! canonical datasets and the replication draws come from identical code.
clear all

capture program drop nns_simdata
program define nns_simdata
    // Generates one simulated example dataset in memory (does not save).
    // design(education|water|clustered) picks the data-generating process.
    // seed(#) makes the draw reproducible; omit to continue the RNG stream.
    syntax, DESign(string) [SEED(integer -1)]

    if !inlist("`design'", "education", "water", "clustered") {
        di as error "design() must be education, water, or clustered"
        exit 198
    }
    if `seed' >= 0 set seed `seed'
    clear

    if "`design'" == "education" {
        // Smooth S-curve with diminishing returns: flat at low spending,
        // rising through the middle, plateauing near the top.
        // True curve: 40 + 45 / (1 + exp(-0.8*(spend - 11))), noise sd = 4.
        set obs 250
        gen district_id = _n
        gen spend_pupil = 5 + 15 * runiform()
        gen true_staar = 40 + 45 / (1 + exp(-0.8 * (spend_pupil - 11)))
        gen staar_pass = true_staar + rnormal(0, 4)
        replace staar_pass = 100 if staar_pass > 100
        replace staar_pass = 0 if staar_pass < 0
        label variable spend_pupil "Per-Pupil Operational Spending ($ Thousands)"
        label variable staar_pass "STAAR Math Passing Rate (%)"
        label variable true_staar "True passing-rate curve (simulation)"
        label data "Texas School District Operational Funding and Achievement (Simulated)"
    }
    else if "`design'" == "water" {
        // Threshold behavior: loss rate is flat below a pipe age of 35 years,
        // then rises with a linear plus quadratic term.  Noise sd = 2.5.
        // The kink at age 35 is the "true change point" referenced in the
        // manuscript's knot-location figure.
        set obs 300
        gen system_id = _n
        gen pipe_age = 80 * runiform()
        gen true_loss = 8 if pipe_age <= 35
        replace true_loss = 8 + 0.4 * (pipe_age - 35) + 0.01 * (pipe_age - 35)^2 ///
            if pipe_age > 35
        gen water_loss = true_loss + rnormal(0, 2.5)
        replace water_loss = 100 if water_loss > 100
        replace water_loss = 0 if water_loss < 0
        label variable pipe_age "Average Water Pipe Age (Years)"
        label variable water_loss "Annual Water Loss Rate (%)"
        label variable true_loss "True loss-rate curve (simulation)"
        label data "Texas Municipal Water Infrastructure and Water Loss (Simulated)"
    }
    else {
        // Clustered support with sharp local features: observations are
        // concentrated in three regions (centers 22, 36, 58) holding a local
        // peak, a sharp dip, and a partial recovery.  This design mirrors the
        // cluster benchmarks in Vinod and Viole (2018), where partition-based
        // fits are expected to have an advantage over a global bandwidth.
        set obs 260
        gen unit_id = _n
        gen mix = runiform()
        gen cluster_score = cond(mix < .45, rnormal(22, 3), ///
            cond(mix < .75, rnormal(36, 2.2), rnormal(58, 4)))
        replace cluster_score = max(5, min(70, cluster_score))
        gen true_index = 18 + .18 * cluster_score ///
            + 10 * exp(-((cluster_score - 24)^2) / 35) ///
            - 16 * exp(-((cluster_score - 35)^2) / 18) ///
            + 8 * exp(-((cluster_score - 58)^2) / 50)
        gen outcome_index = true_index + rnormal(0, 1.2)
        drop mix
        label variable cluster_score "Clustered Index Score"
        label variable outcome_index "Outcome Index"
        label variable true_index "True outcome curve (simulation)"
        label data "Clustered peak/dip benchmark data (Simulated)"
    }
end

// Canonical example datasets used in the manuscript (fixed seeds).
nns_simdata, design(education) seed(2036)
save "examples/texas_education.dta", replace

nns_simdata, design(water) seed(2037)
save "examples/texas_water.dta", replace

nns_simdata, design(clustered) seed(2038)
save "examples/clustered_peak_dip.dta", replace

di "Data simulation complete. Files saved in examples/ subfolder."
