/*==============================================================================
  Lewbel (2012) Instrumental Variables Estimation
  ─────────────────────────────────────────────────
  "Using Heteroscedasticity to Identify and Estimate Mismeasured and
   Endogenous Regressor Models"
  Journal of Business & Economic Statistics, 30(1), 67-80.

  This do-file demonstrates the Lewbel approach to constructing instruments
  from heteroscedasticity when traditional external instruments are
  unavailable or weak.

  Key idea: If the first-stage errors are heteroscedastic, valid instruments
  can be constructed as: Z_j = (X_j - mean(X_j)) * e_hat
  where e_hat are residuals from the first-stage OLS and X_j are exogenous
  regressors.

  Author : Adrij Chakraborty
  Created: 2026-02-12
==============================================================================*/

clear all
set more off

*───────────────────────────────────────────────────────────────────────────────
* 0. Install required packages (run once)
*───────────────────────────────────────────────────────────────────────────────

* ivreg2h   : Lewbel's heteroscedasticity-based IV estimation
* ivreg2    : Extended IV/2SLS/GMM estimation
* ranktest  : Rank test for IV relevance (dependency of ivreg2)

capture which ivreg2h
if _rc {
    ssc install ivreg2h, replace
}
capture which ivreg2
if _rc {
    ssc install ivreg2, replace
}
capture which ranktest
if _rc {
    ssc install ranktest, replace
}

*===============================================================================
* PART A: USING THE ivreg2h COMMAND (RECOMMENDED)
*===============================================================================

*───────────────────────────────────────────────────────────────────────────────
* 1. Load example data
*───────────────────────────────────────────────────────────────────────────────

sysuse auto, clear

* Suppose we want to estimate:
*   price = b0 + b1*mpg + b2*weight + b3*foreign + e
*
* where mpg is endogenous (e.g., correlated with unobserved car quality)
* and we lack conventional external instruments.

describe price mpg weight foreign
summarize price mpg weight foreign

*───────────────────────────────────────────────────────────────────────────────
* 2. OLS baseline (inconsistent if mpg is endogenous)
*───────────────────────────────────────────────────────────────────────────────

regress price mpg weight foreign
estimates store OLS

*───────────────────────────────────────────────────────────────────────────────
* 3. Lewbel IV using ivreg2h — generated instruments only
*───────────────────────────────────────────────────────────────────────────────

* Syntax: ivreg2h depvar exog_vars (endog_vars = [external_IVs]), [options]
*
* When no external instruments are provided, ivreg2h constructs instruments
* solely from heteroscedasticity in the first-stage residuals.

ivreg2h price weight foreign (mpg = )
estimates store Lewbel_only

*───────────────────────────────────────────────────────────────────────────────
* 4. Lewbel IV combined with external instruments
*───────────────────────────────────────────────────────────────────────────────

* If you do have some external instruments, Lewbel instruments can supplement
* them to improve efficiency or to provide overidentification tests.
* Example: suppose 'headroom' and 'trunk' are candidate external instruments.

ivreg2h price weight foreign (mpg = headroom trunk)
estimates store Lewbel_plus_ext

*───────────────────────────────────────────────────────────────────────────────
* 5. Lewbel IV with GMM (efficient under heteroscedasticity)
*───────────────────────────────────────────────────────────────────────────────

ivreg2h price weight foreign (mpg = ), gmm2s
estimates store Lewbel_GMM

*───────────────────────────────────────────────────────────────────────────────
* 6. Compare estimates
*───────────────────────────────────────────────────────────────────────────────

estimates table OLS Lewbel_only Lewbel_plus_ext Lewbel_GMM, ///
    b(%9.3f) se(%9.3f) stats(N r2 j jp)


*===============================================================================
* PART B: MANUAL IMPLEMENTATION (FOR UNDERSTANDING)
*===============================================================================

* This section constructs Lewbel instruments step by step so you can see
* exactly what is happening under the hood.

sysuse auto, clear

*───────────────────────────────────────────────────────────────────────────────
* Step 1: Run first-stage OLS of the endogenous variable on exogenous vars
*───────────────────────────────────────────────────────────────────────────────

regress mpg weight foreign
predict double ehat, residuals

*───────────────────────────────────────────────────────────────────────────────
* Step 2: Test for heteroscedasticity (CRITICAL assumption)
*───────────────────────────────────────────────────────────────────────────────

* The Lewbel approach requires heteroscedasticity in the first-stage errors.
* Without it, the generated instruments have no identifying power.

* Breusch-Pagan / Cook-Weisberg test
quietly regress mpg weight foreign
estat hettest

* White's test (more general)
quietly regress mpg weight foreign
estat imtest, white

display _newline
display as text "────────────────────────────────────────────────────────────"
display as text "NOTE: Reject H0 of homoscedasticity => Lewbel approach is"
display as text "      applicable. If you CANNOT reject, the instruments"
display as text "      will be weak and the method should not be used."
display as text "────────────────────────────────────────────────────────────"

*───────────────────────────────────────────────────────────────────────────────
* Step 3: Construct Lewbel instruments
*───────────────────────────────────────────────────────────────────────────────

* Z_j = (X_j - mean(X_j)) * ehat   for each exogenous regressor X_j

* De-mean the exogenous regressors
foreach var in weight foreign {
    quietly summarize `var'
    generate double `var'_dm = `var' - r(mean)
}

* Construct the generated instruments
generate double z_weight  = weight_dm  * ehat
generate double z_foreign = foreign_dm * ehat

label variable z_weight  "Lewbel IV: (weight - mean) * ehat"
label variable z_foreign "Lewbel IV: (foreign - mean) * ehat"

summarize z_weight z_foreign

*───────────────────────────────────────────────────────────────────────────────
* Step 4: Use the generated instruments in IV estimation
*───────────────────────────────────────────────────────────────────────────────

* 2SLS with generated instruments
ivregress 2sls price weight foreign (mpg = z_weight z_foreign)
estimates store Manual_2SLS

* GMM with generated instruments (robust to heteroscedasticity)
ivregress gmm price weight foreign (mpg = z_weight z_foreign), ///
    wmatrix(robust)
estimates store Manual_GMM

*───────────────────────────────────────────────────────────────────────────────
* Step 5: Post-estimation diagnostics
*───────────────────────────────────────────────────────────────────────────────

* --- First-stage relevance ---
* The generated instruments must be correlated with the endogenous variable
* after partialling out the included exogenous regressors.

regress mpg z_weight z_foreign weight foreign
testparm z_weight z_foreign

display _newline
display as text "────────────────────────────────────────────────────────────"
display as text "First-stage F-statistic for generated instruments."
display as text "Rule of thumb: F > 10 suggests instruments are not weak."
display as text "────────────────────────────────────────────────────────────"

* --- Overidentification test (Hansen J) ---
* Valid only when the model is overidentified (more instruments than
* endogenous regressors). Here we have 2 instruments for 1 endogenous var.

quietly ivregress gmm price weight foreign (mpg = z_weight z_foreign), ///
    wmatrix(robust)
estat overid

display _newline
display as text "────────────────────────────────────────────────────────────"
display as text "Hansen J test of overidentifying restrictions."
display as text "Fail to reject H0 => instruments are plausibly valid."
display as text "────────────────────────────────────────────────────────────"


*===============================================================================
* PART C: MULTIPLE ENDOGENOUS REGRESSORS
*===============================================================================

* The Lewbel approach extends to multiple endogenous variables.
* Each endogenous variable gets its own set of generated instruments.

sysuse auto, clear

* Suppose both mpg and headroom are endogenous:

ivreg2h price weight foreign (mpg headroom = )
estimates store Multi_endog

* With GMM:
ivreg2h price weight foreign (mpg headroom = ), gmm2s
estimates store Multi_endog_GMM

estimates table Multi_endog Multi_endog_GMM, b(%9.3f) se(%9.3f) stats(N j jp)


*===============================================================================
* PART D: PANEL DATA APPLICATION
*===============================================================================

* For panel data, the Lewbel approach can be applied after removing fixed
* effects via within-transformation or first-differencing.

* Example with simulated panel data:

clear
set seed 12345
set obs 500

* Panel structure: 100 individuals, 5 periods
generate id = ceil(_n / 5)
bysort id: generate t = _n
xtset id t

* Generate data with endogeneity
generate double alpha_i = rnormal(0, 1) if t == 1
bysort id (t): replace alpha_i = alpha_i[1]

generate double x1 = rnormal(2, 1)                    // exogenous
generate double x2 = rnormal(0, 1) + 0.5 * alpha_i    // endogenous
generate double u  = rnormal(0, 1) * (1 + 0.5 * abs(x1))  // heteroscedastic

* DGP: y = 1 + 2*x1 + 3*x2 + alpha_i + u
generate double y = 1 + 2*x1 + 3*x2 + alpha_i + u

* --- Within-transform to remove fixed effects ---
foreach var in y x1 x2 {
    bysort id: egen double mean_`var' = mean(`var')
    generate double `var'_w = `var' - mean_`var'
}

* --- Apply Lewbel to within-transformed data ---

* First stage on within-transformed data
regress x2_w x1_w
predict double ehat_w, residuals

* Test heteroscedasticity
estat hettest

* Construct instruments
summarize x1_w
generate double z_x1_w = (x1_w - r(mean)) * ehat_w

* IV estimation on within-transformed data (no constant after FE removal)
ivregress 2sls y_w x1_w (x2_w = z_x1_w), noconstant
display "True coefficients: b1 = 2, b2 = 3"

* GMM version
ivregress gmm y_w x1_w (x2_w = z_x1_w), wmatrix(robust) noconstant
display "True coefficients: b1 = 2, b2 = 3"


*===============================================================================
* SUMMARY
*===============================================================================

display _newline(2)
display as text "================================================================"
display as text "  LEWBEL (2012) IV ESTIMATION — SUMMARY"
display as text "================================================================"
display as text ""
display as text "  When to use:"
display as text "    - Endogenous regressors but no valid external instruments"
display as text "    - External instruments exist but are weak (supplement them)"
display as text "    - Need overidentification tests in just-identified models"
display as text ""
display as text "  Key requirements:"
display as text "    1. Heteroscedasticity in the first-stage equation"
display as text "       (test with Breusch-Pagan or White's test)"
display as text "    2. E[X * e^2] =/= 0 (regressors correlated with error"
display as text "       variance, guaranteed by heteroscedasticity)"
display as text "    3. E[X * e * u] = 0 (product of first- and second-stage"
display as text "       errors uncorrelated with regressors — satisfied if"
display as text "       errors are independent conditional on X)"
display as text ""
display as text "  Quick syntax (ivreg2h):"
display as text "    ivreg2h y x1 x2 (endog = [ext_ivs]), [gmm2s]"
display as text ""
display as text "  Reference:"
display as text "    Lewbel, A. (2012). Using Heteroscedasticity to Identify"
display as text "    and Estimate Mismeasured and Endogenous Regressor Models."
display as text "    JBES, 30(1), 67-80."
display as text "================================================================"

