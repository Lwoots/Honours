#Co-kriging
#Honours 2018
#14 May
##Technical Note: Co-kriging with the gstat package of the R environment for statistical computing, D G Rossiter

#---------------------------------------------------------------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, RColorBrewer, sp, lattice, munsell, gstat)

plots <- read.csv("geo_ref_plots.csv")
sieved <- read.csv("Soil_texture_Final.csv", sep = ";") 
ph <- read.csv("pH_for_analysis.csv", sep = ";")
raw_els <- read.csv("Elsenburg_data_for_analysis.csv", sep = ";")



raw <- left_join(sieved[-c(2:10, 13:14)], ph[-c(2, 6:7)], by = "plot")
raw$munsell[raw$munsell == "10YR 7/3"] <- "10YR 7/4"
raw$munsell[raw$munsell == "7.5YR 6/3"] <- "7.5YR 6/4"


my_data <- left_join(raw, raw_els, by = "plot")


break_down <- mnsl2hvc(my_data$munsell)

my_dat <- data.frame(my_data, break_down)
my_dat$hue <- as.factor(my_dat$hue)
full_dat <- my_dat %>% mutate(my_hue = if_else(hue == "5YR", 5, if_else(hue == "7.5YR", 7.5, if_else(hue == "10YR", 10, 0))))
full_dat <- left_join(plots, full_dat, by = "plot")
grid_data <- full_dat %>% filter(type == "grid")
test_data <- full_dat %>% filter(type == "random")

glimpse(my_dat)
rm(my_data, break_down, raw, sieved, ph, my_dat, raw_els)

#Predicting Mg concentrations by co kriging with soil colour

#The target variable must be normally distributed

hist(grid_data$Mg)
hist(log(grid_data$Mg))
hist(sqrt(grid_data$Mg))
hist(1/grid_data$Mg)

shapiro.test(grid_data$Mg)
shapiro.test(log(grid_data$Mg))
shapiro.test(sqrt(grid_data$Mg)) #Square root transformation is the best

grid_data <- mutate(grid_data, MgSqrt = sqrt(Mg))

#The co variates also need to be normal

hist(grid_data$my_hue)
hist(sqrt(grid_data$my_hue))
shapiro.test(grid_data$my_hue)
shapiro.test(sqrt(grid_data$my_hue)) #Not normal

hist(grid_data$value)
hist(sqrt(grid_data$value))
shapiro.test(grid_data$value)
shapiro.test(sqrt(grid_data$value)) #Not normal

hist(grid_data$chroma)
hist(sqrt(grid_data$chroma))
hist(log(grid_data$chroma))
shapiro.test(grid_data$chroma)
shapiro.test(sqrt(grid_data$chroma)) #Not normal

#That's not great. But let's continue regardless

#Convert datasets into spatial form

coordinates(full_dat) <- c("lon", "lat")
coordinates(grid_data) <- c("lon", "lat")
coordinates(test_data) <- c("lon", "lat")

ggplot(as.data.frame(full_dat), aes(lon, lat)) +
  geom_point(size = sqrt(full_dat$Mg))

#Fit a variogram with no covariates

Mgsq_var <- variog