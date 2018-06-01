rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, RColorBrewer, sp, lattice, munsell, gstat, automap, cowplot)

source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")
rm(species_df)


my_soil_df <- my_soil_df %>% mutate(MgSqrt = sqrt(Mg))
my_grid <- my_soil_df %>% filter(type == "grid")

hist(my_soil_df$ph_kcl)
cor(my_soil_df$ph_kcl, my_soil_df$value, use = "complete")
plot(my_soil_df$ph_kcl ~ my_soil_df$value)

coordinates(my_soil_df) <- c("lon", "lat")
coordinates(my_grid) <- c("lon", "lat")


var_ph <- variogram(ph_kcl ~ 1, my_grid) #Calculate variances between all plots
plot(var_ph)

vgm_ph <- vgm(1.6, "Sph", 0.0003, 0.4) #Create model

fit_ph <- fit.variogram(var_ph, vgm_ph) #fit the model
plot(var_ph, model = fit_ph)

#Do the same for value

var_val <- variogram(value ~ 1, my_grid) #Calculate variances between all plots
plot(var_val)

vgm_val <- vgm(0.28, "Sph", 0.0004, 0.12)

fit_val <- fit.variogram(var_val, vgm_val)
plot(var_val, model = fit_val)

fit_ph$range[2];fit_val$range[2] #Apparently the ranges must be similar - not sure they are though


#Create a gstat object

g <- gstat(NULL, id = "ph", form = ph_kcl ~ 1, data = my_grid) 
g <- gstat(g, id = "value", form = value ~ 1, data = my_soil_df)

var_cross <- variogram(g) #cross variogram
plot(var_cross, pl = T) #ranges don't look too similar

g <- gstat(g, id = "ph", model = fit_ph, fill.all = T) #Fill all gives all variograms the same model, which seems a bit flawed.

g <- fit.lmc(var_cross, g) #fit linear model of co-regionalization
plot(variogram(g), model = g$model) #This looks terrible, but not sure how to fix it


krig <- predict.gstat(g, newdata = my_soil_df)
??predict.gstat

