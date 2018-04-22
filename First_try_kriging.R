#First attempt at kriging
#Honours 2018
#22 April
#Resources: https://rpubs.com/nabilabd/118172, http://rstudio-pubs-static.s3.amazonaws.com/80464_9156596afb2e4dcda53e3650a68df82a.html

#------------------------------------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, RColorBrewer, sp, gstat)

ph <- read.csv("pH_for_analysis.csv", sep = ";")
raw <- read.csv("Species_abundance_for_analysis.csv", sep = ";")

#Extract the names of all the plots
plots <- data.frame(raw[,2])
names(plots) <- "plot"
rm(raw)

#Add pseudo coordinates and site names to the full dataset
full_grid <-  plots %>% mutate(y = as.vector(rep(sort(c(rep(1,11), rep(2,11), rep(3,11), rep(4,11), rep(5,11), rep(6,11), rep(7,11), rep(8,11), rep(9,11), rep(10,11), rep(11,11)), decreasing = T), 3)),
                             x = as.vector(rep(rep(seq(1, 11, 1), 11), 3)), site = c(rep("site1", 121), rep("site2", 121), rep("site3",121)))

glimpse(full_grid)
plot(full_grid$x, full_grid$y)

#now get the smaller grid by joining ph to the full dataset

measured <- left_join(ph, full_grid, by = "plot")
glimpse(measured)

#Extact the full grid of plots from site1
site1_full_grid <- full_grid %>% filter(site == "site1")

#Extract small grid from site 1

site1_measured <- measured %>% filter(site == "site1", type == "grid")

#