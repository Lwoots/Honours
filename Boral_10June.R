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

mod_fit <- fitted(mod1)
pred <- as.data.frame(mod_fit$out)
plot(jitter(pred$tylecodon_pygmaeus) ~ jitter(all_dat$tylecodon_pygmaeus))

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


all_dat_abn <- na.omit(all_dat_abn)

 
sp <- all_dat_abn %>% select(1:15) %>% as.matrix()
covar <- all_dat_abn %>% select(16:35) %>% as.matrix()

grid_abn_allsp_8Aug <- boral(
  sp,
  X = covar,
  family = "negative.binomial",
  num.lv = 3,
  save.model = T
)

summary(grid_abn_allsp_8Aug)
lvsplot(grid_abn_allsp_8Aug)
coefsplot("acidity", grid_abn_allsp_8Aug)


envcors <- get.enviro.cor(grid_abn_allsp_8Aug)
rescors <- get.residual.cor(grid_abn_allsp_8Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor)
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)

calc.varpart(grid_abn_allsp_8Aug)

plot.boral(grid_abn_allsp_8Aug)

grid_abn_allsp_8Aug$geweke.diag$prop.exceed

plot(envcors$cor, rescors$correlation)
plot(envcors$sig.cor, rescors$sig.correlaton)


save(grid_abn_allsp_8Aug, file = "grid_abn_allsp_8Aug.rda")

#predict random plots

test_dat <- my_soil_df %>% 
  filter(type == "random") %>%
  select(c(13:21, 27:28, 30:31, 34, 37:42)) %>% 
  na.omit() %>% 
  as.matrix()

which(is.na(test_dat))

pred <- predict.boral(grid_abn_allsp_8Aug, newX = covar)

grid_abn_allsp_8Aug$call

hist(pred$linpred)
ep <- exp(pred$linpred)
hist(ep)

sum(species_df$ruschia_burtoniae, na.rm = T)

predict(grid_abn_allsp_8Aug, test_dat, type = "response")


#fitted

fit_abn <- fitted.boral(grid_abn_allsp_8Aug)
summary(fit_abn$out)

hist(fit_abn$out)

plot()

pred <- as.data.frame(fit_abn$out)
plot(pred$c_spissum ~ all_dat_abn$c_spissum)
abline(0, 1)
m1 <- lm(pred$c_spissum ~ all_dat_abn$c_spissum)
summary(m1)

pred <- as.data.frame(fit_mod1$out)
plot(pred$c_spissum ~ all_dat$c_spissum)
abline(0, 1)


#redo with non correlated vars ####

enviro_var <- my_soil_df %>% 
  filter(type == "grid") %>%
  select(-c( 
            type, 
            site, 
            munsell, 
            hue,
            ph_h20,
            Very_fine_sand,
            Very_coarse_sand,
            Coarse_sand,
            Medium_sand,
            Fine_sand,
            dN,
            dC,
            Clay,
            Silt,
            acidity,
            conductivity_ms,
            Mg,
            C_perc,
            corr_dC,
            K)) 

all_dat_abn <- left_join(nb_sp_grid_abn, enviro_var, by = "plot") %>% select(-plot)

all_dat_abn <- na.omit(all_dat_abn)


sp <- all_dat_abn %>% select(1:15) %>% as.matrix()
covar <- all_dat_abn %>% select(16:36) %>% as.matrix()

grid_abn_allsp_no_colin_10Aug <- boral(
  sp,
  X = covar,
  family = "negative.binomial",
  num.lv = 2,
  save.model = T
)

save(grid_abn_allsp_no_colin_10Aug, file = "grid_abn_allsp_no_colin_10Aug.rda")

summary(grid_abn_allsp_no_colin_10Aug)


plot.boral(grid_abn_allsp_no_colin_10Aug)

lvsplot(grid_abn_allsp_no_colin_10Aug)
coefsplot("ph_kcl", grid_abn_allsp_no_colin_10Aug)


envcors <- get.enviro.cor(grid_abn_allsp_no_colin_10Aug)
rescors <- get.residual.cor(grid_abn_allsp_no_colin_10Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor)
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)

mod_fit <- fitted.boral(grid_abn_allsp_no_colin_10Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)
plot(pred$c_spissum ~ all_dat_abn$c_spissum)
abline(0, 1)

plot(pred$galenia_fruticosa ~ all_dat_abn$galenia_fruticosa)

which(pred$oophytum == max(pred$oophytum))

grid_abn_allsp_no_colin_10Aug$geweke.diag



#Smaller subset ####

sp <- all_dat_abn %>% select(c(2:8, 10:11, 13:14)) %>% as.matrix()
covar <- all_dat_abn %>% select(16:19, 23:36) %>% as.matrix()


subset_10Aug <- boral(
  sp,
  X = covar,
  family = "negative.binomial",
  num.lv = 2,
  save.model = T
)

save(subset_10Aug, file = "subset_10Aug.rda")
plot.boral(subset_10Aug)

summary(subset_10Aug)

mod_fit <- fitted.boral(subset_10Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)
plot(pred$c_spissum ~ all_dat_abn$c_spissum)
plot(pred$a_delaetii ~ all_dat_abn$a_delaetii)
plot(pred$oophytum[-13] ~ all_dat_abn$oophytum[-13])
plot(pred$c_subfenestratum ~ all_dat_abn$c_subfenestratum)
plot(pred$a_fissum ~ all_dat_abn$a_fissum)

abline(0,1)



#Conservative try ####

species_df[is.na(species_df)] <- 0

spdat <- species_df %>% left_join(my_soil_df, by = "plot") %>% filter(!is.na(type))

spdat <- na.omit(spdat)
sp <- spdat %>% select(ruschia_burtoniae, c_spissum, a_fissum, a_delaetii, c_subfenestratum, dicrocaulon, mesemb_1, drosanthemum_diversifolium, crassula_muscosa, tylecodon_pygmaeus, oophytum, a_framesii) %>% as.matrix()

covar <- spdat %>% select(acidity, Sand, drainage, Olsen, Na, C_N_ratio, Ca) %>% as.matrix()


simple_mod <- boral(
  sp,
  X = covar,
  family = "negative.binomial",
  num.lv = 2,
  save.model = T
)

plot.boral(simple_mod)

mod_fit <- fitted.boral(simple_mod, est = "mean")
pred <- as.data.frame(mod_fit$out)
plot(pred$drosanthemum_diversifolium ~ spdat$drosanthemum_diversifolium)

summary(pred)
abline(0,1)

mod_fit$ordinal.probs

save(simple_mod, file =  "simple_mod.rda")
