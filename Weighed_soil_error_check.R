setwd("/Users/larawootton/Documents/Honours/Data")
rm(list = ls())

### Checking the raw weighed soil texture data for errors ###

dat <- read.csv("Soil_texture_data21Mar_explore.csv", sep = ";")

summary(dat)


dotchart(dat$total_pan)
max(dat$total_pan)
hist(dat$total_pan)
boxplot(dat$total_pan)

dotchart(dat$sieve_weight)
max(dat$sieve_weight)
hist(dat$sieve_weight)
boxplot(dat$sieve_weight)

dotchart(dat$soil_mass)
max(dat$soil_mass)
hist(dat$soil_mass)
boxplot(dat$soil_mass)

dotchart(dat$less2mm)
max(dat$less2mm)
hist(dat$less2mm)
boxplot(dat$less2mm)


dotchart(dat$between2_5mm)
max(dat$between2_5mm)
hist(dat$between2_5mm)
boxplot(dat$between2_5mm)

dotchart(dat$total_less2mm)
max(dat$total_less2mm)
hist(dat$total_less2mm)
boxplot(dat$total_less2mm)

dotchart(dat$less1mm)
max(dat$less1mm)
hist(dat$less1mm)
boxplot(dat$less1mm)

dotchart(dat$between2_1)
max(dat$between2_1)
hist(dat$between2_1)
boxplot(dat$between2_1)





