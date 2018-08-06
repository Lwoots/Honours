#Joint species distribution modelling
#Getting to grips with boral
#Honours 2018
#10 June

#---------------------------------------------------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr,visdat, ggplot2, RColorBrewer, boral, rjags, mvabund)


source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

species_df <- species_df %>% select(-c(Q_cover, X, site, transect, lampranthus))


species_df <- species_df %>% select(1:28)
sp_pres <- species_df[,-c(1:3)] %>% apply(2, function(x) ifelse(is.na(x), 0, 1))


#Covariates
my_soil_df <- my_soil_df %>% filter(!is.na(type))

covar <- as.matrix(my_soil_df[,c(13:21, 27:28, 30:31, 34, 37)])

mod1 <- boral(sp_pres, X = covar, family = "binomial")
summary(mod1)
lvsplot(mod1)




mod2 <- boral(sp_pres, num.lv = 3, family = "binomial")
lvsplot(mod2)
coefsplot(mod2)

#Spider tutorial


data("spider")

y <- spider$abund

fit.lvmp <- boral(y = y, family = "poisson", num.lv = 2, row.eff = "fixed")
X <- scale(spider$x)
fit.Xnb <- boral(y = y, X = X, family = "negative.binomial", num.lv = 2)



species_df <- species_df %>% select(-c(Q_cover, X, site, transect, lampranthus))
b <- species_df[4:16]
glimpse(b)
mod3 <- boral(y = b, family = "poisson", num.lv = 2, row.eff = "fixed")
plot(mod3, ask = F, mfrow = c(2,2))
lvsplot(mod3)


dat <- left_join(species_df, my_soil_df, by = c("plot", "lat", "lon", "site"))
dat <- dat %>% filter(!is.na(type))
dat <- dat %>% select(-c(Q_cover, X, site, transect, lampranthus))

b <- dat[,4:22]

v <- scale(dat[,c(58:64, 73, 75)])

mod4 <- boral(y = b, X = v, family = "negative.binomial", num.lv = 3, save.model = T)
plot(mod4)
lvsplot(mod4)

envcors <- get.enviro.cor(mod4)
rescors <- get.residual.cor(mod4)


library(corrplot)
corrplot(envcors$cor, type = "lower", diag = F, title = "Correlations due
         to covariates", mar = c(3,0.5,2,1), tl.srt = 45)
corrplot(rescors$cor, type = "lower", diag = F, title = "Residual correlati
ons", mar = c(3,0.5,2,1), tl.srt = 45)


# 22 July ####

sp <- species_df %>% select(ruschia_burtoniae, c_spissum, a_delaetii, c_subfenestratum, dicrocaulon, drosanthemum_diversifolium, plot) 

vars <- my_soil_df %>% select(acidity, Mg, Clay, K, C_N_ratio, conductivity_ms, lon, lat, drainage, plot, type) %>% na.omit

testdat <- left_join(sp, vars, by = "plot") %>% filter(!is.na(type)) %>% select(-type)

y <- as.matrix(testdat[,1:6])
X <- as.matrix(testdat[,8:16])

my_boral_abn <- boral(y, 
                      X, 
                      family = "negative.binomial", 
                      num.lv = 2,
                      save.model = T
                      )
glimpse(my_boral_abn)

summary(my_boral_abn)
plot(my_boral_abn)

lvsplot(my_boral_abn)


rescorr <- get.residual.cor(my_boral_abn)
corrplot(rescorr$correlation, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)

corrplot(rescorr$sig.correlaton, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)


envcorr <- get.enviro.cor(my_boral_abn)

corrplot(envcorr$cor)


#Occurrence

y_occ <- ifelse(y>0,1,0)

my_boral_occ <- boral(y_occ, 
                      X, 
                      family = "binomial", 
                      num.lv = 2,
                      save.model = T
)


summary(my_boral_occ)
plot(my_boral_occ)

lvsplot(my_boral_occ)


rescorr <- get.residual.cor(my_boral_occ)
corrplot(rescorr$correlation, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)

corrplot(rescorr$sig.correlaton, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)


envcorr <- get.enviro.cor(my_boral_occ)

corrplot(envcorr$cor)
corrplot(envcorr$sig.cor)


my_boral_abn$geweke.diag$prop.exceed
my_boral_abn