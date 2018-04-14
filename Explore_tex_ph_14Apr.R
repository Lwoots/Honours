#Honours 2018
#Exploratory analysis 14 April
#Descriptive stats for texture and ph, conductivity

#----------------------------------------------------#
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, visdat)

plots <- read.csv("plots.csv")
sieved <- read.csv("Soil_texture_Final.csv", sep = ";") 
malvern <- read.csv("Raw_Malvern_data_for_analysis_31Mar.csv", sep = ";")

which(duplicated(malvern$Sample) == TRUE) #Find the duplicated samples in the malvern data

malvern <- malvern[-c(156, 153, 152, 65, 81, 95, 99, 103),]



?unique
duplicated

setdiff(plots$plot, sieved$plot)
setdiff(plots$plot, malvern$Sample)
setdiff(malvern$Sample, sieved$plot)
