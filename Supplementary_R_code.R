#######################################################################################
#######################################################################################
### Time-stratified case-crossover studies in environmental epidemiology:           ###
### a tutorial.                                                                     ###
### (2nd version submitted to IJE Educational Corner, 2023/04/01)                   ### 
#######################################################################################
#######################################################################################

# Remove objects.
rm(list = ls())

# Install packages.
install.packages("dplyr", "splines", "gnm", "Epi")

# Load packages.
library(dplyr); library(splines); library(gnm); library(Epi) 

#######################################################################################
# Data management.
#######################################################################################

list.cities <- list()
list.cities[[1]]<-read.csv('valencia.csv')
list.cities[[2]]<-read.csv('london.csv')
cities <- c('valencia', 'london')

for (city in 1:length(list.cities))
  {
  
  # Load raw datset.
  data <- list.cities[[city]]
  
  # Generate lags for temperature, humidity, and PM10.
  for (num in 0:3)
  {
    tmean <- paste("tmean_l", num, sep="")
    assign(tmean, dplyr::lag(data$tmean, num))
    rh <- paste("rh_l", num, sep="")
    assign(rh, dplyr::lag(data$rh, num))
    pm10 <- paste("pm10_l", num, sep="")
    assign(pm10, dplyr::lag(data$pm10, num))
  }
  
  # Generate splines for temperature adjustment.
  data$l03tmean <- (tmean_l0+tmean_l1+tmean_l2+tmean_l3)/4
  data$ns.tmean <- ns(data$l03tmean, df=6)
  
  # Generate splines for humidity adjustment.
  data$l03rh <- (rh_l0+rh_l1+rh_l2+rh_l3)/4
  data$ns.rh <- ns(data$l03rh, df=3)
  
  # Generate averaged lag exposure for PM10 sale for a 10 ug/m3 increase.
  data$l01pm10 <- (pm10_l0+pm10_l1)/2
  data$l01pm10 <- data$l01pm10/10

	# Save R data file.
	save(data, file = paste0(cities[city],".RData"))
}

#######################################################################################
# Example 1. Exposure-outcome association and environmental time-varying confounders 
# adjustment using the Valencia dataset.
#######################################################################################

# Load Valencia dataset.
load("valencia.RData")

# Generate time-stratified strata.
data$month <- as.factor(data$month)
data$year  <- as.factor(data$year)
data$dow   <- as.factor(data$dow)
data$stratum <- with(data, as.factor(year:month:dow))
data$ind <- tapply(data$all, data$stratum, sum)[data$stratum]

# Fit fixed-effects conditional quasi-Poisson regression.
model.cc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)

# Get Relative Risk for PM10.
Epi::ci.exp(model.cc, subset="l01pm10") 

#######################################################################################
# Adjustment of subpopulation time-invariant covariates.
#######################################################################################

# Reshape to long format.
long.all <- reshape2::melt(dplyr::select(data, date, allage1, allage2),
                           id.vars="date", variable.name="age", value.name="allage")
long <- merge(long.all, dplyr::select(data, -c(allage1, allage2)), by="date", all.x=T)

# Fit fixed-effects conditional quasi-Poisson regression adjusted by age.
model.cc.adj.age <- gnm(allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
                        data=long, family=quasipoisson, eliminate=stratum)
Epi::ci.exp(model.cc.adj.age, subset="l01pm10") 

# Generate age-time-stratified strata.
long$stratum4 <- with(long, as.factor(year:month:dow:age))
long$ind4 <- tapply(long$all, long$stratum4, sum)[long$stratum4]

# Fit conditional quasi-Poisson with age-time-stratified strata.
model.cc.str4 <- gnm(allage ~ ns.tmean + ns.rh + l01pm10, 
                     data=long, family=quasipoisson, subset=ind4>0, eliminate=stratum4)
Epi::ci.exp(model.cc.str4, subset="l01pm10") 

#######################################################################################
# Investigation of effect modification by age.
#######################################################################################

# Stratified analysis by age. 
model.cc.age1 <- gnm(allage1 ~ ns.tmean + ns.rh + l01pm10, 
                     data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.age1, subset="l01pm10") 

model.cc.age2 <- gnm(allage2 ~ ns.tmean + ns.rh + l01pm10, 
                     data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.age2, subset="l01pm10") 

# Interaction analysis by age. 
model.cc.age.int <- gnm(allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
                        data=long, family=quasipoisson, subset=ind4>0, eliminate=stratum4)
Epi::ci.exp(model.cc.age.int, subset="l01pm10") 

# Likelihood ratio-test for effect for interaction.
anova(model.cc.str4, model.cc.age.int, test="LRT")

#######################################################################################
# Example 2. Multi-location study.
#######################################################################################

# Load Valencia R data file.
load("valencia.RData")
data.vlc <- dplyr::select(data, -allage1,-allage2)

# Load London R data file.
load("london.RData")
data.ldn <- data

# Append datasets.
dat.2city <- rbind(data.vlc, data.ldn)

# Generate space-time-stratified strata.
dat.2city$city  <- as.factor(dat.2city$city)
dat.2city$month <- as.factor(dat.2city$month)
dat.2city$year  <- as.factor(dat.2city$year)
dat.2city$dow   <- as.factor(dat.2city$dow)
dat.2city$stratum <- with(dat.2city, as.factor(city:year:month:dow))
dat.2city$ind <- tapply(dat.2city$all, dat.2city$stratum, sum)[dat.2city$stratum]

# Fit conditional Poisson with space-time-stratified strata.
model.cc.adj.city <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                         data=dat.2city, family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.adj.city, subset="l01pm10") 

# Stratified analysis by city.
model.cc.vlc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data=subset(dat.2city,city=="Valencia"), 
                    family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.vlc, subset="l01pm10") 

model.cc.ldn <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data=subset(dat.2city,city=="London"), 
                    family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.ldn, subset="l01pm10") 

# Interaction analysis by location.
model.cc.city.int <- gnm(all ~ ns.tmean + ns.rh + city + l01pm10:city, 
                        data=dat.2city, family=quasipoisson, subset=ind>0, eliminate=stratum)
Epi::ci.exp(model.cc.city.int, subset="l01pm10") 

# Likelihood ratio test for effect modification.
anova(model.cc.adj.city, model.cc.city.int, test="LRT")

#######################################################################################
#######################################################################################
###                           End of script file                                    ###
#######################################################################################
#######################################################################################
