setwd("/Users/larawootton/Documents/Honours/Data")
rm(list = ls())

### Checking the raw weighed soil texture data for errors ###

dat <- read.csv("Soil_texture_data21Mar_explore.csv", sep = ";")

summary(dat)


dotchart(dat$between2_5mm)
max(dat$between2_5mm)
hist(dat$between2_5mm)
boxplot(dat$between2_5mm)

