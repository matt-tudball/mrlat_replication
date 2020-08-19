clear all
set more off
global path "yourpath"
set maxvar 30000

// Format genetic data
/*
import delimited "$path/DATA/bmi_genetic_raw.raw", delimiter(space) 
destring rs*, force replace
drop iid pat mat sex phenotype
rename fid n_eid
save "$path/DATA/bmi_genetic_raw.dta", replace
*/

// Set up
use "path to phenotype data file", clear

// Merge SNP data
merge 1:1 n_eid using "path to genetic data file"
keep if _merge == 3
drop _merge

// Merge genetic data
merge 1:1 n_eid using "path to principal components file"
keep if _merge == 3
drop _merge

// Generate main variables
// Genetic instruments
local i = 1
foreach var of varlist rs* {
	rename `var' z`i'
	local ++i
}

// Outcome controls
gen c_sex = n_31_0_0 // Sex
gen c_centre = n_54_0_0 // Assessment centre
gen c_age = n_21003_0_0 // Age
gen c_alc = n_1558_0_0 // Alcohol intake (frequency)
gen c_smoke = n_20116_0_0 // Smoking status
gen c_towns = n_189_0_0 // Townsend deprivation index
gen c_batch = n_22000_0_0 + 20 // Genotype measurement batch
forvalues i = 1/10 {
	gen c_pc`i' = n_22009_0_`i' // First 10 principal components
	drop if missing(c_pc`i')
}
drop if missing(c_towns)
drop if missing(c_batch)

// Liability
gen liability = n_21001_0_0 // BMI

// Other liability exposures
// gen x_yob = n_34_0_0 // Year of birth
// gen x_msmoke = n_1787_0_0 // Maternal smoking around birth
// gen x_size10 = n_1687_0_0 // Comparative body size age 10
// gen x_hgt10 = n_1697_0_0 // Comparative height age 10
// gen x_nation = n_1647_0_0 // Country of birth
// gen x_mum = n_1835_0_0 // Mother still alive

// Outcome
// gen outcome = n_40007_0_0 // All-cause mortality
// replace outcome = 0 if missing(outcome)
// replace outcome = 1 if outcome > 0 & !missing(outcome)
gen outcome = n_4080_0_0 // Systolic blood pressure

keep z* c_* liability outcome

foreach var of varlist z* {
	count if missing(`var')
	if r(N) > 1000 {
		drop `var'
		continue
	}
	drop if missing(`var')
}

foreach var of varlist c_sex - c_smoke liability outcome {
	drop if `var' < 0
	drop if missing(`var')
}

saveold "$yourpath/example_data.dta", replace version(12)

// Analyse data
gen diagnosis = (liability >= 30) // Obesity

// First stage prune instruments
reg liability z*, robust
foreach var of varlist z* {
	local tstat = _b[`var']/_se[`var']
	if `tstat' > 4 | `tstat' < -4 {
		rename `var' incl_`var' 
	}
}

// Second stage
global controls_y c_sex c_age i.c_alc i.c_centre i.c_smoke c_towns i.c_batch c_pc* 
global controls_g 

reg liability incl_z* z* $controls_g c_sex, robust 

ivregress 2sls outcome (liability = incl_z*) $controls_y, robust
estat firststage

qui sum liability
gen liability_std = (liability - `r(mean)')/`r(sd)'
ivregress 2sls outcome (liability_std = incl_z*) $controls_y, robust
estat firststage

qui probit diagnosis incl_z* z* $controls_g c_sex
predict Ghat, xb

/*
cumul liability, generate(CDFliability)
ivregress 2sls outcome (CDFliability = incl_z*) $controls_y, robust
estat firststage

cumul Ghat, generate(CDFGhat)
ivregress 2sls outcome (CDFGhat = incl_z*) $controls_y, robust
estat firststage
*/

qui sum Ghat
gen Ghat_std = (Ghat - `r(mean)')/`r(sd)'
ivregress 2sls outcome (Ghat_std = incl_z*) $controls_y, robust
estat firststage

ivregress 2sls outcome (diagnosis = incl_z*) $controls_y, robust
estat firststage

capture program drop ivliab
program ivliab, rclass
	syntax, in1(string) out1(name) out2(name) //out3(name)
	qui {
		local theta = sqrt(`in1')
		probit diagnosis incl_z* z* $controls_g c_sex
		predict Ghat, xb
		
		// Lower bound on standard deviation scale
		qui sum Ghat
		gen Ghat_std = (Ghat - `r(mean)')/`r(sd)'
		ivregress 2sls outcome (Ghat_std = incl_z*) $controls_y, robust
		local beta_sdl = _b[Ghat_std]
		
		// Exact on standard deviation scale
		local beta_sde = `beta_sdl'/`theta'
		
		// Lower bound on percentile rank scale
		/*
		cumul Ghat, generate(CDFGhat)
		ivregress 2sls outcome (CDFGhat = incl_z*) $controls_y, robust
		local beta_prl = _b[CDFGhat]
		*/
		
		drop Ghat Ghat_std // CDFGhat
	}
	return local `out1' = `beta_sdl'
	return local `out2' = `beta_sde'
	// return local `out3' = `beta_prl'
end

/*
ivliab, in1(0.0718) out1(beta_sdl) out2(beta_sde)
di "`r(beta_sde)'"
bootstrap real(r(beta_sdl)) real(r(beta_sde)), reps(10): lower, in1(0.0718) ///
	out1(beta_sdl) out2(beta_sde)
*/
foreach var in Ghat_std Ghat diagnosis {
	capture drop `var'
}
local theta = 0.0256 // Variance explained
local nobs = 6 // Number of cut-offs
local i = 0 // Index value
foreach k in 22.5 25 27.5 30 32.5 35 {
	local i = `i'+1
	
	gen diagnosis = (liability >= `k')

	// ivliab, in1(`sigmag') out1(beta_sdl) out2(beta_sde) // out3(beta_prl)
	bootstrap real(r(beta_sdl)) real(r(beta_sde)), seed(1) reps(1000): ///
		ivliab, in1(`theta') out1(beta_sdl) out2(beta_sde)
	local beta_sdl`i' = _b[_bs_1]
	local beta_sde`i' = _b[_bs_2]
	local std_sdl`i' = _se[_bs_1]
	local std_sde`i' = _se[_bs_2]
	
	ivregress 2sls outcome (diagnosis = incl_z*) $controls, robust
	local beta_d`i' = _b[diagnosis]
	local std_d`i' = _se[diagnosis]
	
	// convolution, in1(`sigmag') in2(50) out1(beta_pre)
	// local beta_pre`i' = `r(beta_pre)'
	
	// convolution, in1(0.1) in2(100) out1(betasl)
	// local betasl`i' = `r(betasl)'
	
	local beta_sdldi = round(`beta_sdl`i'',0.01)
	local beta_sdedi = round(`beta_sde`i'',0.01)
	// local beta_prldi = round(`beta_prl`i'',0.01)
	local beta_ddi = round(`beta_d`i'',0.01)
	// local beta_predi = round(`beta_pre`i'',0.01)
	// local betasldi = round(`betasl`i'',0.01)
	
	di "Cut-off: `k', Lower bound: `beta_sdldi', Diagnosis: `beta_ddi', Exact: `beta_sdedi'"
	
	drop diagnosis 
}

preserve
	clear
	capture set obs `nobs'
	gen cutoff = .
	gen beta_sdl = .
	gen beta_sde = .
	// gen beta_prl = .
	gen beta_d = .
	// gen beta_pre = .
	// gen betasu = .
	gen std_sdl = .
	gen std_sde = .
	gen std_d = .
	local j = 0
	foreach k in 22.5 25 27.5 30 32.5 35 {
		local j = `j'+1
		replace cutoff = `k' in `j'
		replace beta_sdl = `beta_sdl`j'' in `j'
		replace beta_sde = `beta_sde`j'' in `j'
		// replace beta_prl = `beta_prl`j'' in `j'
		replace beta_d = `beta_d`j'' in `j'
		// replace beta_pre = `beta_pre`j'' in `j'
		// replace betasu = `betasu`j'' in `j'
		replace std_sdl = `std_sdl`j'' in `j'
		replace std_sde = `std_sde`j'' in `j'
		replace std_d = `std_d`j'' in `j'
	}
	saveold "$yourpath/vary_cutoff.dta", replace version(12)
	
	// Reshape from wide to long
	rename beta_sdl beta1
	rename beta_sde beta2
	rename beta_d beta3
	rename std_sdl std1
	rename std_sde std2
	rename std_d std3
	
	reshape long beta std, i(cutoff) j(group)
	saveold "$yourpath/vary_cutoff_long.dta", replace version(12)
restore
