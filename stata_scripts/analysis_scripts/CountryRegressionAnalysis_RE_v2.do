/*================================================================================

AUTHOR:					PWS
DATE: 					October 2017
VERSION:				STATA/MP 15.0
DO FILE NAME:			CountryRegressionAnalysis

STATUS:					final

DEPENDENCIES:			MergeVars.do

DATASETS USED: 			complete_build.dta - Final build of dataset, ready for analysis

DESCRIPTION OF FILE:	UK country comparison of audit analysis (NHLI poster)

==================================================================================*/

clear
set more off

cd "C:\Users\pstone\OneDrive - Imperial College London\PhD\Objective 2 & 3 - CPRD NCAP Primary Care Audit Replication\CPRD GOLD Analysis"

/* Create log file */
capture log close
log using analysis_logs/CountryRegressionAnalysis_RE_v2, smcl replace

local data_dir "D:\CPRD NCAP Primary Care Audit Replication\CPRD GOLD\Analysis"

use "`data_dir'/builds/complete_build", clear


//covariates
gen byte agebands = 1 if age >= 35
replace agebands = 2 if age >= 45
replace agebands = 3 if age >= 55
replace agebands = 4 if age >= 65
replace agebands = 5 if age >= 75
replace agebands = 6 if age >= 85
replace agebands = . if age == .
label define agebands 1 "35-44" 2 "45-54" 3 "55-64" 4 "65-74" 5 "75-84" 6 "85+"
label values agebands agebands
order agebands, after(age)

keep if gender == 1 | gender == 2    //include men and women only

tab agebands, missing
tab gender, missing

tab anxiety_bin, missing
tab asthma_bin, missing
tab bronchiectasis_bin, missing
tab coronary_heart_disease_bin, missing
tab depression_bin, missing
tab diabetes_bin, missing
tab heart_failure_bin, missing
tab pain, missing
tab hypertension_bin, missing
tab lung_cancer_bin, missing
tab osteoporosis_bin, missing
tab psychosis_bin, missing
tab stroke_bin, missing

tab mrc_score, missing
tab currentsmokstatus, missing


//exposure
label list country
tab country, missing

gen country2 = 1 if country == 3     //new variable so that Wales is at top
replace country2 = 2 if country == 1
replace country2 = 3 if country == 2
replace country2 = 4 if country == 4
label define country2 1 "Wales" 2 "England" 3 "Scotland" 4 "Northern Ireland"
label values country2 country2

tab country country2


//outcomes
gen byte fev1fvc = lowfev1fvcratio
replace fev1fvc = 0 if fev1fvc == 2

//= EXTRA SPIROMETRY DETAIL ====================================================
tab country2 lowfev1fvcratio, row chi   //339m codes for diagnoses in past 2 yrs
tab country2 anyfev1fvcratio, row chi   //any spirometry codes for diagnoses in past 2 yrs
//==============================================================================

gen byte mrcpastyr = mrc_pastyr
replace mrcpastyr = 1 if mrc_pastyr >= 1

gen byte asksmokpastyr = smokstatus_pastyr
replace asksmokpastyr = 1 if asksmokpastyr >= 1

gen byte mrc3plus_pr3yrs_nox = 0 if mrc_score >= 3 & mrc_score <= 5
replace mrc3plus_pr3yrs_nox = 1 if mrc_score >= 3 & mrc_score <= 5 ///
								& pulmonary_rehab != . ///
								& pulmonary_rehab > end-(3*365.25)

tab fev1fvc                 //Restricted to diagnoses in the last 2 years
//sensitivity: tab anyfev1fvcratio         //Restricted to diagnoses in the last 2 years
tab cxr_6mo                 //Restricted to diagnoses in the last 2 years
tab mrcpastyr, missing
tab asksmokpastyr, missing
tab flu_vaccine, missing
tab behavanddrug            //Restricted to smokers in the last 2 years
tab mrc3plus_pr3yrs         //Restricted to those without exception report for PR
//sensitivity:
tab mrc3plus_pr3yrs_nox, missing



//i) Diagnosis
tab fev1fvc country2, col chi
melogit fev1fvc i.country2 || pracid:, base or
melogit fev1fvc i.country2 i.agebands i.gender || pracid:, base or
//melogit fev1fvc i.country2 i.agebands i.gender mrc_score currentsmokstatus || pracid:, base or
melogit fev1fvc i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
/*melogit fev1fvc i.country2 i.agebands i.gender mrc_score currentsmokstatus ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or*/
estimates store spirometry
estat icc

//tab anyfev1fvcratio country2, col chi
//melogit anyfev1fvcratio i.country2 || pracid:, base or

tab cxr_6mo country2, col chi
melogit cxr_6mo i.country2 || pracid:, base or
melogit cxr_6mo i.country2 i.agebands i.gender || pracid:, base or
melogit cxr_6mo i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store cxr
estat icc


//ii) Assessment
tab mrc_pastyr country, missing col chi   //look at difference in distribution of scores
tab mrc_score country, missing col chi

tab mrcpastyr country2, col chi
melogit mrcpastyr i.country2 || pracid:, base or
melogit mrcpastyr i.country2 i.agebands i.gender || pracid:, base or
melogit mrcpastyr i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store mrc_grade
estat icc

tab asksmokpastyr country2, col chi
melogit asksmokpastyr i.country2 || pracid:, base or
melogit asksmokpastyr i.country2 i.agebands i.gender || pracid:, base or
melogit asksmokpastyr i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store smoking
estat icc


//iii) Treatment
tab flu_vaccine country2, col chi
melogit flu_vaccine i.country2 || pracid:, base or
melogit flu_vaccine i.country2 i.agebands i.gender || pracid:, base or
melogit flu_vaccine i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store fluvax
estat icc

tab behavanddrug country2, col chi
melogit behavanddrug i.country2 || pracid:, base or
melogit behavanddrug i.country2 i.agebands i.gender || pracid:, base or
melogit behavanddrug i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store smokcess
estat icc

tab mrc3plus_pr3yrs country2, col chi
melogit mrc3plus_pr3yrs i.country2 || pracid:, base or
melogit mrc3plus_pr3yrs i.country2 i.agebands i.gender || pracid:, base or
melogit mrc3plus_pr3yrs i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estimates store pr_ref
estat icc

tab mrc3plus_pr3yrs_nox country2, col chi
melogit mrc3plus_pr3yrs_nox i.country2 || pracid:, base or
melogit mrc3plus_pr3yrs_nox i.country2 i.agebands i.gender || pracid:, base or
melogit mrc3plus_pr3yrs_nox i.country2 i.agebands i.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin || pracid:, base or
estat icc


coefplot spirometry cxr mrc_grade smoking fluvax smokcess pr_ref ///
		 , graphregion(color(white)) bgcolor(white) ///
		 xlab(0.1 0.5(0.5)3, grid glcolor(gs15)) ylab(, glcolor(gs15)) ///
		 xtitle("Odds ratio") ///
		 eform xscale(log range(0.02 3)) ///
		 coeflabels(, wrap(12)) ciopts(recast(rcap)) ///
		 xline(1, lcolor(black) lwidth(thin) lpattern(dash)) ///
		 drop(_cons *.agebands *.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin) ///
		 legend(order(2 "Confirmation of airway obstruction" 4 "Chest X-ray" 6 "Record of MRC grade in the past year" 8 "Record of smoking status in the past year" 10 "Receipt of seasonal influenza immunisation" 12 "Smoking cessation treatment" 14 "Referred to pulmonary rehabilitation"))
graph save graphs/all_outcomes, replace


coefplot (spirometry, color(blue) ciopts(lc(blue) recast(rcap))) ///
		 (cxr, color(orange) ciopts(lc(orange) recast(rcap))) ///
		 , graphregion(color(white)) bgcolor(white) ///
		 xlab(0.1 0.5(0.5)3, grid glcolor(gs15)) ylab(, glcolor(gs15)) ///
		 xtitle("Odds ratio") ///
		 eform xscale(log range(0.02 3)) ///
		 coeflabels(, wrap(12)) ///
		 xline(1, lcolor(black) lwidth(thin) lpattern(dash)) ///
		 drop(_cons *.agebands *.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin) ///
		 title("Diagnosis", color(black)) ///
		 legend(region(lcolor(white)) col(1) ///
		 order(2 "Confirmation of airway obstruction" 4 "Chest X-ray" 6 ""))
graph save graphs/spirometry_cxr, replace
		 
coefplot (mrc_grade, color(gray) ciopts(lc(gray) recast(rcap))) ///
		 (smoking, color(gold) ciopts(lc(gold) recast(rcap))) ///
		 , graphregion(color(white)) bgcolor(white) ///
		 xlab(0.1 0.5(0.5)3, grid glcolor(gs15)) ylab(, glcolor(gs15)) ///
		 xtitle("Odds ratio") ///
		 eform xscale(log range(0.02 3)) yscale(off) ///
		 coeflabels(, wrap(12)) ///
		 xline(1, lcolor(black) lwidth(thin) lpattern(dash)) ///
		 drop(_cons *.agebands *.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin) ///
		 title("Assessment", color(black)) ///
		 legend(region(lcolor(white)) col(1) ///
		 order(2 "Record of MRC grade in the past year" 4 "Record of smoking status in the past year" 6 ""))
graph save graphs/mrc_smoking, replace

coefplot (fluvax, color(midblue) ciopts(lc(midblue) recast(rcap))) ///
		 (smokcess, color(green) ciopts(lc(green) recast(rcap))) ///
		 (pr_ref, color(dknavy) ciopts(lc(dknavy) recast(rcap))) ///
		 , graphregion(color(white)) bgcolor(white) ///
		 xlabel(0.1 0.5(0.5)3, grid glcolor(gs15)) ylabel(, glcolor(gs15)) ///
		 xtitle("Odds ratio") ///
		 eform xscale(log range(0.02 3)) yscale(off) ///
		 coeflabels(, wrap(12)) ///
		 xline(1, lcolor(black) lwidth(thin) lpattern(dash)) ///
		 drop(_cons *.agebands *.gender ///
		 anxiety_bin asthma_bin bronchiectasis_bin coronary_heart_disease_bin ///
		 depression_bin diabetes_bin heart_failure_bin pain hypertension_bin ///
		 lung_cancer_bin osteoporosis_bin psychosis_bin stroke_bin) ///
		 title("High-value care", color(black)) ///
		 legend(region(lcolor(white)) col(1) ///
		 order(2 "Receipt of seasonal influenza immunisation" 4 "Smoking cessation treatment" 6 "Referred to pulmonary rehabilitation"))
graph save graphs/fluvax_smokcess_pr, replace


graph combine graphs/spirometry_cxr.gph /// 
			  graphs/mrc_smoking.gph ///
			  graphs/fluvax_smokcess_pr.gph ///
			  , graphregion(color(white)) scale(0.95) cols(3) xsize(30) ysize(10)
graph save graphs/combined_three, replace
graph export graphs/combined_three.svg, replace
graph export graphs/combined_three.pdf, replace
graph export graphs/combined_three.eps, replace


log close
translate analysis_logs/CountryRegressionAnalysis_RE_v2.smcl analysis_logs/CountryRegressionAnalysis_RE_v2.log