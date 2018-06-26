rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")
source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, dismo, gbm)

glimpse(my_soil_df)

#Select important environmental variables
nb_soil_var <- my_soil_df %>% select(plot:lat, type:percent_over1, value:Very_coarse_sand, conductivity_ms, ph_kcl, N_perc, corr_dN:C_perc, corr_dC:Q_cover)
glimpse(nb_soil_var)

#Select important species
nb_sp <- species_df %>% select(1,6:10, 12, 16:17, 20, 22, 24,28:33)
nb_sp[is.na(nb_sp)] <- 0

#Create presence-absence matrix

prez <- ifelse(nb_sp[,2:length(nb_sp)] > 0, 1, 0)
prez_df <- data.frame(plot = nb_sp[,1], prez)
my_prez <- left_join(prez_df, nb_soil_var, by = "plot")
my_prez <- my_prez %>% filter(type == "grid" | type == "random") %>% select(-type)

#Create abundance df
my_data <- left_join(nb_sp, nb_soil_var, by = "plot")
my_data <- my_data %>% filter(type == "grid" | type == "random")
my_data <- my_data %>% select(-type)


#Set BRT parameters
brt_var <- 19:51
response <- 7
tree.com <- 5
learn <- 0.001




brt_results <- gbm.step(my_prez,
                       gbm.x = brt_var,
                       gbm.y = response,
                       plot.main = TRUE,
                       family = "bernoulli",
                       step.size = 50,
                       tree.complexity = tree.com,
                       learning.rate = learn,
                       max.trees=10000,
                       n.folds = 10,
                       bag.fraction = 0.5
                       )


(pseudo_r2 <- 1-(brt_results$cv.statistics$deviance.mean/brt_results$self.statistics$mean.null))

summary(brt_results)
gbm.plot(brt_results)



brt.simp <- gbm.simplify(brt_results, n.drops = "auto")
summary(brt.simp)
