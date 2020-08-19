clear all
set more off
global path "yourpath"
set maxvar 30000

// Format genetic data
/*
import delimited "path to genetic .raw file (females)", delimiter(space) 
destring rs*, force replace
drop iid pat mat sex phenotype
rename fid n_eid
save "save genetic data", replace
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

//Program to combine the fields 0 and 1 for each variable
cap prog drop merge_var
prog def merge_var
	replace `1'_0_0=`1'_1_0 if (`1'_0_0==.|`1'_0_0<0 )& (`1'_1_0>0 & `1'_1_0!=.)
end

// Generate main variables
foreach var of varlist rs* {
	count if missing(`var')
	if r(N) > 2000 {
		drop `var'
		continue
	}
	drop if missing(`var')
}

// Keep females
keep if n_31_0_0 == 0

// Outcome controls
// gen c_centre = n_54_0_0 // Assessment centre
gen c_age = n_21003_0_0 // Age
gen c_batch = n_22000_0_0 + 20 // Genotype measurement batch
/*
forvalues i = 1/10 {
	gen c_pc`i' = n_22009_0_`i' // First 10 principal components
	drop if missing(c_pc`i')
}
*/

// Dichotomous exposure
gen d_size = n_1687_0_0 // Comparative body size age 10

// Outcome (doctor-diagnosed diabetes)
/*
merge_var n_2443
merge_var n_2976
gen out_diabetes=(n_2443_0_0==1) if n_2443_0_0>=0 
replace out_diabete=. if n_2976_0_0<=21 & n_2976_0_0!=.
*/ 

keep rs* c_* d_*

foreach var of varlist c_* d* {
	drop if `var' < 0
	drop if missing(`var')
}

saveold "$yourpath/richardson_data_females.dta", replace version(12)

// Ordered probit regression
oprobit d_size rs* c_age i.c_batch
foreach var of varlist rs* {
	local `var'_beta = _b[`var']
	local `var'_se = _se[`var']
	qui sum `var'
	local `var'_eaf = r(mean)/2
	local `var'_var = r(Var)
}

keep rs*
foreach var of varlist rs* {
	drop if !missing(`var')
}

set obs 4

foreach var of varlist rs* {
	replace `var' = ``var'_beta' in 1
	replace `var' = ``var'_se' in 2
	replace `var' = ``var'_eaf' in 3
	replace `var' = ``var'_var' in 4
}

xpose, varname clear
rename _varname SNP
rename v1 beta
rename v2 se
rename v3 eaf
rename v4 var

split SNP, p("_")
drop SNP
rename SNP1 SNP
rename SNP2 effect_allele
replace effect_allele = upper(effect_allele)
gen other_allele = "T" if effect_allele == "A"
replace other_allele = "A" if effect_allele == "T"
replace other_allele = "C" if effect_allele == "G"
replace other_allele = "G" if effect_allele == "C"

saveold "$yourpath/richardson_twosample_females.dta", replace version(12)
outsheet using "$yourpath/richardson_twosample_females.csv", replace comma 
