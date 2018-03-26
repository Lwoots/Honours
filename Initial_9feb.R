rm(list = ls())

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, visdat)

## Exploring the data 9 Feb ##

setwd("/Users/larawootton/Documents/Honours/Data")

kners <- read.csv("Knersvlakte_20Jan_explore.csv", sep = ";")
kners$site <- as.factor(kners$site)
head(kners)
summary(kners)

vis_dat(kners[,3:51])


str(kners$galenia_fruticosa)
summary(kners$galenia_fruticosa)
summary(kners$Q_cover)

kners$galenia_fruticosa[kners$galenia_fruticosa == "cf 1"] <- 1
kners$galenia_fruticosa <- as.integer(kners$galenia_fruticosa)
str(kners$monilaria_moniliformis)


plot(kners$site, kners$a_delaetii)
plot(kners$site, kners$c_spissum)
plot(kners$site, kners$a_fissum)

site1 <- kners[kners$site == "1",]
site2 <- kners[kners$site == "2",]
site3 <- kners[kners$site == "3",]

site1$transect <- factor(site1$transect)
site2$transect <- factor(site2$transect)
site3$transect <- factor(site3$transect)

par(mfrow = c(1,3))
plot(site1$transect, site1$a_delaetii)
plot(site2$transect, site2$a_delaetii)
plot(site3$transect, site3$a_delaetii)

levels(site1$transect) <- c("AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK")

summary(site3)

pairs(kners[4:10])
