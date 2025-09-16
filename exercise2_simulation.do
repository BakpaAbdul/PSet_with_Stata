/**************************************************************************************************
  Problem set 2 — Randomized Experiment Simulation (Stata .do)
  Author: Abdul Baari Bakpa
  Purpose: simulating randomized experiments to see how the plain diff-in-means (the ATE estimator) actually behaves in practice—its bias, variance, and CI coverage—under controlled data-generating processes.
           Part A (Q6):   Fixed finite population; re-randomize treatment 800 times.
           Part B (Q7):   Redraw a new population each repetition.
**************************************************************************************************/

version 19.5
clear all
set more off

// --------------------------- USER PARAMETERS --------------------------------
local N        = 1000          // population size
local Ntreat   = 500           // number treated each assignment
local reps     = 800           // Monte Carlo repetitions
local seedA    = 12345         // seed for Part A (fixed pop, re-assign)
local seedB    = 54321         // seed for Part B (redraw population each rep)
// ---------------------------------------------------------------------------

// --------------------------- HELPER: ASSERTIONS -----------------------------
program define _assert_vars, rclass
    // Ensure required variables exist
    syntax varlist(min=1)
    foreach v of varlist `varlist' {
        capture confirm variable `v'
        if _rc {
            di as error "Required variable `v' not found."
            exit 111
        }
    }
end

// --------------------------- PART A: FIXED POP ------------------------------
program drop _all
program define ex2_partA, rclass
    version 18.0
    syntax, N(integer) Ntreat(integer) reps(integer) seed(integer)

    // 1) Build finite population once
    clear
    set seed `seed'
    set obs `N'
    gen double y0 = invnormal(runiform())              // U ~ N(0,1)
    gen double v  = invnormal(runiform())              // V ~ N(0,1), indep. of U
    gen double y1 = 0.5*y0 + 0.5*v + 0.2              // potential outcome under treatment
    gen double te = y1 - y0

    quietly summarize te, meanonly
    scalar ATE_N = r(mean)
    di as text "Part A — Finite-population ATE_N = " as result %9.6f ATE_N

    // 2) Re-randomize treatment; collect estimates
    tempfile resultsA
    tempname pfA
    postfile `pfA' double(b v cover) using "`resultsA'", replace

    forvalues r = 1/`reps' {
        quietly {
            gen double rand = runiform()
            sort rand
            gen byte D  = _n <= `Ntreat'                // exactly Ntreat treated
            gen double Y = (1-D)*y0 + D*y1

            regress Y D, vce(robust)
            scalar b  = _b[D]
            scalar vv = _se[D]^2

            scalar lo = b - 1.96*sqrt(vv)
            scalar hi = b + 1.96*sqrt(vv)
            scalar cover = (ATE_N >= lo) & (ATE_N <= hi)

            post `pfA' (b) (vv) (cover)
            drop rand D Y
        }
    }
    postclose `pfA'

    use "`resultsA'", clear
    quietly summarize b, meanonly
    scalar mean_bA = r(mean)
    scalar var_bA  = r(Var)
    quietly summarize v, meanonly
    scalar mean_vA = r(mean)
    quietly summarize cover, meanonly
    scalar covA = r(mean)

    di as text "Part A — mean(b) vs ATE_N: " as result %9.6f mean_bA as text " vs " as result %9.6f ATE_N
    di as text "Part A — Var(b): " as result %9.6f var_bA as text "   mean robust Vhat: " as result %9.6f mean_vA
    di as text "Part A — 95% CI coverage for ATE_N: " as result %6.3f covA
    // keep results in memory for user; return scalars
    return scalar ATE_N = ATE_N
    return scalar mean_b = mean_bA
    return scalar var_b  = var_bA
    return scalar mean_v = mean_vA
    return scalar cover  = covA
end

// --------------------------- PART B: REDRAW POP -----------------------------
program define ex2_partB, rclass
    version 18.0
    syntax, N(integer) Ntreat(integer) reps(integer) seed(integer)

    // Target super-population constant effect component is 0.2
    local tau0 = 0.2

    tempfile resultsB
    tempname pfB
    postfile `pfB' double(b v cover02) using "`resultsB'", replace

    set seed `seed'
    forvalues r = 1/`reps' {
        quietly {
            clear
            set obs `N'
            gen double y0 = invnormal(runiform())
            gen double v  = invnormal(runiform())
            gen double y1 = 0.5*y0 + 0.5*v + `tau0'

            gen double rand = runiform()
            sort rand
            gen byte D  = _n <= `Ntreat'
            gen double Y = (1-D)*y0 + D*y1

            regress Y D, vce(robust)
            scalar b  = _b[D]
            scalar vv = _se[D]^2
            scalar lo = b - 1.96*sqrt(vv)
            scalar hi = b + 1.96*sqrt(vv)
            scalar cover02 = (`tau0' >= lo) & (`tau0' <= hi)

            post `pfB' (b) (vv) (cover02)
        }
    }
    postclose `pfB'

    use "`resultsB'", clear
    quietly summarize b, meanonly
    scalar mean_bB = r(mean)
    scalar var_bB  = r(Var)
    quietly summarize v, meanonly
    scalar mean_vB = r(mean)
    quietly summarize cover02, meanonly
    scalar covB = r(mean)

    di as text "Part B — mean(b) vs 0.2: " as result %9.6f mean_bB as text " vs 0.2"
    di as text "Part B — Var(b): " as result %9.6f var_bB as text "   mean robust Vhat: " as result %9.6f mean_vB
    di as text "Part B — 95% CI coverage for 0.2: " as result %6.3f covB

    return scalar mean_b = mean_bB
    return scalar var_b  = var_bB
    return scalar mean_v = mean_vB
    return scalar cover  = covB
end

// ============================ RUN THE PARTS ==================================
quietly capture log close _all
// Uncomment to log: log using "exercise2_simulation.smcl", replace

// Part A
ex2_partA, N(`N') Ntreat(`Ntreat') reps(`reps') seed(`seedA')

// Part B
ex2_partB, N(`N') Ntreat(`Ntreat') reps(`reps') seed(`seedB')

// End of file
