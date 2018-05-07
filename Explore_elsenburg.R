#Exploratory analysis of elsenburg soil data from the knersvlakte fieldtrip January
#Honours 2018
#3 May

#
#Set up------------------------------------------------------------------------
#

rm(list = ls())

setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, visdat, ggplot2, RColorBrewer, cowplot, corrplot)

raw_els <- read.csv("Elsenburg_data_for_analysis.csv", sep = ";")
glimpse(raw_els)

plots <- read.csv("plots.csv")


setdiff(plots$plot, raw_els$plot)
setdiff(raw_els$plot, plots$plot)

grid_els <- raw_els %>% arrange(plot)

site3 <- grid_els[7:42,] %>% 
  mutate(y = sort(as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6))), decreasing = T),
         x = as.vector(rep(seq(1, 11, 2), 6)))

site1 <- grid_els[c(1:6, 43:72),] %>% 
  mutate(y = sort(as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6))), decreasing = T),
         x = as.vector(rep(seq(1, 11, 2), 6)))

temp <- matrix(rep(NA, length(grid_els)), nrow = 1)
colnames(temp) <- colnames(grid_els)
site2 <- rbind(grid_els[c(73:90, 93:96),], temp, grid_els[c(97, 100:111),]) %>% 
  mutate(y = sort(
    as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6))), decreasing = T),
    x = as.vector(rep(seq(1, 11, 2), 6)))
rm(temp)



#exploratory analysis####
vis_dat(grid_els)
pairs(grid_els)
c <- cor(grid_els[,2:8])
corrplot.mixed(c, lower.col = "black")


#acidity####

par(mfrow = c(1,3))
hist(grid_els$acidity, xlab = "Exchangeable acidity (cmol/kg)", main = "")
boxplot(grid_els$acidity, ylab = "Exchangeable acidity (cmol/kg)")
points(mean(grid_els$acidity), pch = 19, col = "blue", cex = 1.2)
dotchart(grid_els$acidity, lcolor = "white", ylab = "Order of data", xlab = "Exchangeable acidity (cmol/kg)")
