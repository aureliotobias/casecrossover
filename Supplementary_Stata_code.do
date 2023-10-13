********************************************************************************
********************************************************************************
*** Time-stratified case-crossover studies in environmental epidemiology:    ***
*** a tutorial.                                                              ***
*** (2nd version submitted to IJE Educational Corner, 2023/04/01)            ***
********************************************************************************
********************************************************************************


* Install Stata user-written commands. 
net install gensplines , from("https://www.pclambert.net/downloads/gensplines") force
net install xtodp , from("https://raw.githubusercontent.com/aureliotobias/conditionalpoisson/master/") force

* Note: gensplines only works in Stata 17.

********************************************************************************
* Data management.
********************************************************************************

foreach city in valencia london {
	* Load raw datset.
	import delimited "`city'.csv" , clear
	
	* Encoding date.
	generate date2 = date(date, "DMY")
	format date2 %td
	tsset date2 , daily
	
	* Generate lags for temperature, humidity, and PM10.
	foreach num of numlist 0/3 {
		generate tmean_l`num' = L`num'.tmean
		generate rh_l`num' = L`num'.rh
		generate pm10_l`num' = L`num'.pm10
	}
	
	* Generate splines for temperature adjustment.
	generate l03tmean = (tmean_l0+tmean_l1+tmean_l2+tmean_l3)/4
	gensplines l03tmean , type(ns) df(6) gen(nstmean_)

	* Generate splines for humidity adjustment.
	generate l03rh = (rh_l0+rh_l1+rh_l2+rh_l3)/4
	gensplines l03rh , type(ns) df(3) gen(nsrh_)
	
	* Generate averaged lag exposure for PM10 sale for a 10 ug/m3 increase.
	generate l01pm10 = (pm10_l0+pm10_l1)/2
	replace l01pm10 = l01pm10/10
	
	* Save Stata dataset.
	save "`city'" , replace
}

********************************************************************************
* Example 1. Exposure-outcome association and time-varying confounders 
* adjustment using the Valencia dataset.
********************************************************************************

* Load Valencia dataset.
use "valencia" , replace

* Generate time-stratified strata.
egen stratum = group(year month dow)
xtset stratum

* Fit fixed-effects conditional Poisson regression and correct by overdispersion.
xtpoisson all nstmean_* nsrh_* l01pm10 , fe
xtodp

* Get Relative Risk for PM10.
lincom l01pm10 , eform

********************************************************************************
* Adjustment of subpopulation time-invariant covariates.
********************************************************************************

* Reshape to long format.
reshape long allage , i(date) j(age)

* Fit conditional Poisson (and correct by overdispersion) adjusted by age.
xtpoisson allage nstmean_* nsrh_* i.age l01pm10 , fe
estimates store A1
xtodp
lincom l01pm10 , eform

* Generate age-time-stratified strata.
egen stratum4 = group(year month dow age)
xtset stratum4

* Fit conditional Poisson (and correct by overdispersion) with age-time-stratified strata.
xtpoisson allage nstmean_* nsrh_* i.age l01pm10 , fe
estimates store A
xtodp
lincom l01pm10 , eform

********************************************************************************
* Investigation of effect modification by age.
********************************************************************************

* Stratified analysis by age.
xtset stratum
foreach num of numlist 1 2 {
	xtpoisson allage nstmean_* nsrh_* l01pm10 if age==`num' , fe
	xtodp if age==`num'
	lincom l01pm10 , eform
}

* Interaction analysis by age.
xtset stratum4
xtpoisson allage nstmean_* nsrh_* c.l01pm10##i.age , fe
estimates store B
xtodp
lincom l01pm10 , eform                     
lincom (l01pm10 + 2.age#l01pm10) , eform   

* Likelihood ratio-test for effect for interaction.
lrtest A B

********************************************************************************
* Example 2. Multi-location study.
********************************************************************************

* Load Valencia dataset.
use "valencia" , replace

* Append London data set.
append using "london"

* Generate space-time-stratified strata.
encode city , generate(citynum)
egen stratum = group(year month dow citynum)
xtset stratum

* Fit conditional Poisson (and correct by overdispersion) with space-time-stratified strata.
xtpoisson all nstmean_* nsrh_* c.l01pm10 , fe
estimates store A
xtodp
lincom l01pm10 , eform

* Stratified analysis by location.
foreach name in Valencia London  {
	xtpoisson all nstmean_* nsrh_* l01pm10 if city=="`name'" , fe
	xtodp if city=="`name'"
	lincom l01pm10 , eform
}

* Interaction analysis by location.
xtpoisson all nstmean_* nsrh_* c.l01pm10##i.citynum , fe
estimates store B
xtodp
lincom l01pm10 , eform                         
lincom (l01pm10 + 2.citynum#l01pm10) , eform   

* Likelihood ratio test for effect modification.
lrtest A B

********************************************************************************
********************************************************************************
***                            End of do file                                ***
********************************************************************************
********************************************************************************
