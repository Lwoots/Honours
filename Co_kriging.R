#Co-kriging
#Honours 2018
#14 May
##Technical Note: Co-kriging with the gstat package of the R environment for statistical computing, D G Rossiter

#---------------------------------------------------------------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, RColorBrewer, sp, lattice, munsell, gstat, automap, cowplot)

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
test_data <- mutate(as.data.frame(test_data), MgSqrt = sqrt(Mg))

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
full_dat <- mutate(as.data.frame(full_dat), MgSqrt = sqrt(Mg))
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

Mgsq_var_cloud <- variogram(MgSqrt ~ 1, data = grid_data, cloud = T)
plot(Mgsq_var_cloud)

Mgsq_var <- variogram(MgSqrt ~ 1, data = grid_data)
plot(Mgsq_var, pl = T)

#fit model to variogram

Mgsq_vgm <- vgm(0.65, "Sph", range =  0.00037, 0.18) #This is total thumb suck. Just played around until it fitted
plot(Mgsq_var, pl = T, model = Mgsq_vgm)

#This isn't looking hopeful, so let's auto fit.

auto_Mg <- autofitVariogram(MgSqrt ~ 1, grid_data)

summary(auto_Mg)
plot(auto_Mg)

autokrige_Mg <- autoKrige(MgSqrt ~ 1, grid_data, full_dat)
plot(autokrige_Mg)
summary(autokrige_Mg)

p <- ggplot(as.data.frame(full_dat), aes(lon, lat, colour = MgSqrt)) +
  geom_point(size = 5) +
  theme_bw()

pred_krige_mg <- as.data.frame(autokrige_Mg[1])
p2 <- ggplot(pred_krige_mg, aes(krige_output.lon, krige_output.lat, colour = krige_output.var1.pred)) +
  geom_point(size = 5) +
  theme_bw()

plot(pred_krige_mg)

plot_grid(p, p2)

full_dat <- as.data.frame(full_dat)
colnames(pred_krige_mg)[1] <- c("lon")
pred <- left_join(full_dat, pred_krige_mg, by = "lon")

test <- pred %>% filter(type == "random")

plot(test$krige_output.var1.pred ~ test$MgSqrt)
mod <- lm(test$krige_output.var1.pred ~ test$MgSqrt)
summary(mod) #r squared = 0.59
abline(mod)

#Now for co kriging


coordinates(full_dat) <- c("lon", "lat")
coordinates(grid_data) <- c("lon", "lat")
coordinates(test_data) <- c("lon", "lat")


co_auto_Mg <- autofitVariogram(MgSqrt ~ chroma, grid_data)

summary(auto_Mg)
plot(auto_Mg)

co_autokrige_Mg <- autoKrige(MgSqrt ~ chroma, grid_data, full_dat)
plot(co_autokrige_Mg)
summary(co_autokrige_Mg)


p <- ggplot(as.data.frame(full_dat), aes(lon, lat, colour = MgSqrt)) +
  geom_point(size = 5) +
  theme_bw()

pred_krige_mg <- as.data.frame(co_autokrige_Mg[1])
p2 <- ggplot(pred_krige_mg, aes(krige_output.lon, krige_output.lat, colour = krige_output.var1.pred)) +
  geom_point(size = 5) +
  theme_bw()

plot(pred_krige_mg)

plot_grid(p, p2)

#three = 0.2053




full_dat <- as.data.frame(full_dat)
colnames(pred_krige_mg)[1] <- c("lon")
pred <- left_join(full_dat, pred_krige_mg, by = "lon")

test <- pred %>% filter(type == "random")

plot(test$krige_output.var1.pred ~ test$MgSqrt)
mod <- lm(test$krige_output.var1.pred ~ test$MgSqrt)
summary(mod) #r squared = 0.59
abline(mod)


plot(pred$MgSqrt ~ pred$chroma)




auto_Mg <- autofitVariogram(Mg ~ 1, grid_data)

summary(auto_Mg)
plot(auto_Mg)

autokrige_Mg <- autoKrige(Mg ~ 1, grid_data, full_dat)
plot(autokrige_Mg)
summary(autokrige_Mg)

p <- ggplot(as.data.frame(full_dat), aes(lon, lat, colour = MgSqrt)) +
  geom_point(size = 5) +
  theme_bw()

pred_krige_mg <- as.data.frame(autokrige_Mg[1])
p2 <- ggplot(pred_krige_mg, aes(krige_output.lon, krige_output.lat, colour = krige_output.var1.pred)) +
  geom_point(size = 5) +
  theme_bw()

plot(pred_krige_mg)

plot_grid(p, p2)

full_dat <- as.data.frame(full_dat)
colnames(pred_krige_mg)[1] <- c("lon")
pred <- left_join(full_dat, pred_krige_mg, by = "lon")

test <- pred %>% filter(type == "random")

plot(test$krige_output.var1.pred ~ test$MgSqrt)
mod <- lm(test$krige_output.var1.pred ~ test$Mg)
summary(mod) #r squared = 0.59
abline(mod)