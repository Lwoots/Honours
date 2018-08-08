#Joint species distribution modelling
#Getting to grips with boral
#Honours 2018
#10 June

#---------------------------------------------------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr,visdat, ggplot2, RColorBrewer, boral, rjags, mvabund, corrplot)


source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

nb_sp <- species_df %>% left_join(my_soil_df, by = "plot") %>% filter(!is.na(type)) %>%  select(plot, ruschia_burtoniae, c_spissum, a_delaetii, oophytum, a_fissum, dicrocaulon, drosanthemum_diversifolium, a_framesii, mesemb_1, c_subfenestratum, co_calculus, galenia_fruticosa, crassula_muscosa, tylecodon_pygmaeus, cephalophyllum_staminodiosum) 

enviro_var <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  select(c(1, 13:21, 27:28, 30:31, 34, 37:42))


sp_occ <- nb_sp[,2:16] %>% apply(2, function(x) ifelse(is.na(x), 0, 1))

all_dat <- cbind(sp_occ, enviro_var)
all_dat <- na.omit(all_dat) %>% select(-plot)

#Full model ####

#Covariates
covar <- all_dat %>% 
  select(c(16:35)) %>% 
  as.matrix()

sp <- as.matrix(all_dat[,1:15])

mod1 <- boral(sp, 
              X = covar,
              family = "binomial",
              num.lv = 3,
              save.model = T)

summary(mod1)
lvsplot(mod1)
coefsplot("Q_cover", mod1)


envcors <- get.enviro.cor(mod1)
rescors <- get.residual.cor(mod1)

corrplot(envcors$cor)
corrplot(envcors$sig.cor)
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)
fitted.boral(mod1)
var <- calc.varpart(mod1)
plot(var$varpart.X)

fitted(mod1)
plot.boral(mod1)


save(mod1, file = "mod1_6Aug.rda")
load("mod1_6Aug.rda")

#Model without sandy species ####

covar <- all_dat %>% 
  select(c(16:35)) %>% 
  as.matrix()

sp <- as.matrix(all_dat[,c(1:8, 10:14)])

less_sandsp_7Aug <- boral(sp, 
              X = covar,
              family = "binomial",
              num.lv = 3,
              save.model = T)
summary(less_sandsp_7Aug)
lvsplot(less_sandsp_7Aug)
coefsplot("Q_cover", less_sandsp_7Aug)


envcors <- get.enviro.cor(less_sandsp_7Aug)
rescors <- get.residual.cor(less_sandsp_7Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor)
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)

calc.varpart(less_sandsp_7Aug)

save(less_sandsp_7Aug, file = "less_sandsp_7Aug.rda")

#Quartz species only ####

sp <- as.matrix(all_dat[,c(2:8, 10:11, 13:14)])

quartz_7Aug <- boral(
  sp,
  X = covar,
  family = "binomial",
  num.lv = 3,
  save.model = T
)
summary(quartz_7Aug)
lvsplot(quartz_7Aug)
coefsplot("acidity", quartz_7Aug)


envcors <- get.enviro.cor(quartz_7Aug)
rescors <- get.residual.cor(quartz_7Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor)
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)

calc.varpart(quartz_7Aug)

quartz_7Aug$geweke.diag$prop.exceed

quartz_7Aug$prior.control

save(quartz_7Aug, file = "quartz_7Aug.rda")

#Abundance all sp based on grid plots ####

nb_sp_grid_abn <- species_df %>%
  left_join(my_soil_df, by = "plot") %>%
  filter(type == "grid") %>%
  select(
    plot,
    ruschia_burtoniae,
    c_spissum,
    a_delaetii,
    oophytum,
    a_fissum,
    dicrocaulon,
    drosanthemum_diversifolium,
    a_framesii,
    mesemb_1,
    c_subfenestratum,
    co_calculus,
    galenia_fruticosa,
    crassula_muscosa,
    tylecodon_pygmaeus,
    cephalophyllum_staminodiosum
  ) 

nb_sp_grid_abn[is.na(nb_sp_grid_abn)] <- 0 #replace NA with zero

enviro_var <- my_soil_df %>% 
  filter(type == "grid") %>%
  select(c(1, 13:21, 27:28, 30:31, 34, 37:42)) 

all_dat_abn <- left_join(nb_sp_grid_abn, enviro_var, by = "plot") %>% select(-plot)




 


enviro_var <- my_soil_df %>% 
  filter(type == "grid") %>%
  select(c(1, 13:21, 27:28, 30:31, 34, 37:42)) %>% 
  as.matrix()

grid_abn_allsp_8Aug <- boral(
  nb_sp_grid_abn,
  X = enviro_var,
  family = "negative.binomial",
  num.lv = 3,
  save.model = T
)

which(is.na(enviro_var))

#Ordination

ord <- boral(sp,
              family = "binomial",
              num.lv = 2,
              save.model = T)

lvsplot(ord)



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
