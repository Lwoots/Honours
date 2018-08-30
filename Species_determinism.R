#Niche determinism and specificity
#http://www.will.chez-alice.fr/pdf/KleyerJVS2012AppendixS2.pdf p23

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(ggplot2, adehabitatHS, MASS, dplyr, boral, foreach, doParallel, TeachingDemos, vegan, rlist)

#Setup ####
source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")


soil_dat <- my_soil_df %>% 
  select(site, 
         ph_kcl,
         Na,
         C_N_ratio,
         Ca,
         Q_cover,
         elevation,
         percent_over1,
         percent_over2,
         P,
         N_perc,
         aspect,
         Clay,
         corr_dN,
         drainage,
         slope,
         type) %>% 
  filter(!is.na(type)) %>% 
  na.omit() %>% select(-type)

soil_dat[,2:16] <- scale(soil_dat[,2:16])

#Select study species (those that are present in ten or more plots) 
nb_sp <- species_df %>% select(site,
                               plot,
                               ruschia_burtoniae,
                               mesemb_1, 
                               drosanthemum_diversifolium,
                               a_delaetii,
                               a_fissum,
                               a_framesii,
                               c_spissum,
                               cephalophyllum_staminodiosum,
                               dicrocaulon,
                               oophytum) 


names(nb_sp) <- c("site", 
                  "plot",
                  "R_burtoniae", 
                  "R_comptonii", 
                  "D_diversifolium",
                  "A_delaetii",
                  "A_fissum",
                  "A_framesii",
                  "C_spissum",
                  "C_staminodiosum",
                  "Dicrocaulon_sp",
                  "Oophytum_sp")

nb_sp[is.na(nb_sp)] <- 0 #replace NA with 0

sp_occ <- nb_sp[3:12] %>% 
  apply(2, function(x) ifelse(x == 0, 0, 1)) %>% 
  as.data.frame() %>% 
  mutate(site = nb_sp$site) %>% 
  mutate(plot = nb_sp$plot) %>% 
  na.omit()

all_dat <- left_join(sp_occ, my_soil_df %>% select(colnames(soil_dat[2:15]), plot), by = "plot") %>% 
  na.omit()

all_dat[,13:26] <- scale(all_dat[,13:26])

rm(nb_sp, sp_occ, soil_dat)

#Stochastic null model ####

stoc_auc_burt <- list()
stoc_auc_comp <- list()
stoc_auc_div <- list()
stoc_auc_del <- list()
stoc_auc_fis <- list()
stoc_auc_fram <- list()
stoc_auc_spis <- list()
stoc_auc_stam <- list()
stoc_auc_dic <- list()
stoc_auc_ooph <- list()
stoc_sp_num <- list()

#set up parallel
cores <- detectCores()
cl<-makeCluster(cores - 1)
registerDoParallel(cl)

num = 1

repeat{
  set.seed(Sys.time())
  rplots <- sample(150, 100)
  train <- all_dat[rplots, ]
  test <-  all_dat[-c(rplots), ]
  
  #Shuffle species, keeping row and column sums equal to original data.
  rsp <- permatswap(train[,1:10], times = 1, burnin = 5, fixedmar = "both", mtype = "prab")
  
  #Record the number of species used in each plot
  stoc_sp_num[[num]] <- train %>% summarise(R_burtoniae = sum(R_burtoniae),
                                             R_comptonii = sum(R_comptonii),
                                             D_diversifolium = sum(D_diversifolium),
                                             A_delaetii = sum(A_delaetii),
                                             A_fissum = sum(A_fissum),
                                             A_framesii = sum(A_framesii),
                                             C_spissum = sum(C_spissum),
                                             C_staminodiosum = sum(C_staminodiosum),
                                             Dicrocaulon_sp = sum(Dicrocaulon_sp),
                                             Oophytum_sp = sum(Oophytum_sp)
  )
  
  #Skip loop interation
  
  if(stoc_sp_num[[num]][[8]] == 10) {
    next
  }
  
  sp <- as.data.frame(rsp$perm)
  covar <- as.matrix(train[,13:length(train)])
  
  mod <- boral(
    sp,
    covar,
    num.lv = 3,
    family = "binomial",
    #mcmc.control = example_mcmc_control,
    save.model = T
  )
  
  test_covar <-  as.matrix(test[,13:length(test)]) #Note that the test data here is unpermuted
  
  newpred <- predict.boral(mod, 
                           newX = test_covar, 
                           predict.type = "marginal",
                           est = "mean")
  #Calculate auc for each species
  stoc_auc_burt[[num]] <- pROC::roc(test$R_burtoniae, newpred$linpred[,1]) %>% pROC::auc()
  stoc_auc_comp[[num]] <- pROC::roc(test$R_comptonii, newpred$linpred[,2]) %>% pROC::auc()
  stoc_auc_div[[num]] <- pROC::roc(test$D_diversifolium, newpred$linpred[,3]) %>% pROC::auc()
  stoc_auc_del[[num]] <- pROC::roc(test$A_delaetii, newpred$linpred[,4]) %>% pROC::auc()
  stoc_auc_fis[[num]] <- pROC::roc(test$A_fissum, newpred$linpred[,5]) %>% pROC::auc()
  stoc_auc_fram[[num]] <- pROC::roc(test$A_framesii, newpred$linpred[,6]) %>% pROC::auc()
  stoc_auc_spis[[num]] <- pROC::roc(test$C_spissum, newpred$linpred[,7]) %>% pROC::auc()
  stoc_auc_stam[[num]] <- pROC::roc(test$C_staminodiosum, newpred$linpred[,8]) %>% pROC::auc()
  stoc_auc_dic[[num]] <- pROC::roc(test$Dicrocaulon_sp, newpred$linpred[,9]) %>% pROC::auc()
  stoc_auc_ooph[[num]] <- pROC::roc(test$Oophytum_sp, newpred$linpred[,10]) %>% pROC::auc()
  
   
  num <- num + 1
  if(num > 10) {
    break
  }
}

stoc_run1_n2_28Aug <- list(stoc_auc_burt,
                           stoc_auc_comp,
                           stoc_auc_div ,
                           stoc_auc_del ,
                           stoc_auc_fis ,
                           stoc_auc_fram,
                           stoc_auc_spis,
                           stoc_auc_stam,
                           stoc_auc_dic ,
                           stoc_auc_ooph)

stoc_run2_n10_28Aug <- list(stoc_auc_burt,
                           stoc_auc_comp,
                           stoc_auc_div ,
                           stoc_auc_del ,
                           stoc_auc_fis ,
                           stoc_auc_fram,
                           stoc_auc_spis,
                           stoc_auc_stam,
                           stoc_auc_dic ,
                           stoc_auc_ooph)

save(stoc_run1_n2_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/stoc_run1_n2_28Aug.rda")
save(stoc_run2_n10_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/stoc_run2_n10_28Aug.rda")

sp_stoc_run1_n2_28Aug <- sp_num
save(sp_stoc_run1_n2_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/sp_stoc_run1_n2_28Aug.rda")


stoc_R_burt <- c(unlist(stoc_run1_n2_28Aug[[1]]), unlist(stoc_run2_n10_28Aug[[1]]))
stoc_R_comp <- c(unlist(stoc_run1_n2_28Aug[[2]]),unlist(stoc_run2_n10_28Aug[[2]]))
stoc_D_div <- c(unlist(stoc_run1_n2_28Aug[[3]]),unlist(stoc_run2_n10_28Aug[[3]]))
stoc_A_del <- c(unlist(stoc_run1_n2_28Aug[[4]]),unlist(stoc_run2_n10_28Aug[[4]]))
stoc_A_fis <- c(unlist(stoc_run1_n2_28Aug[[5]]),unlist(stoc_run2_n10_28Aug[[5]]))
stoc_A_fra <- c(unlist(stoc_run1_n2_28Aug[[6]]),unlist(stoc_run2_n10_28Aug[[6]]))
stoc_C_spis <- c(unlist(stoc_run1_n2_28Aug[[7]]),unlist(stoc_run2_n10_28Aug[[7]]))
stoc_C_sta <- c(unlist(stoc_run1_n2_28Aug[[8]]),unlist(stoc_run2_n10_28Aug[[8]]))
stoc_Dicro <- c(unlist(stoc_run1_n2_28Aug[[9]]),unlist(stoc_run2_n10_28Aug[[9]]))
stoc_Ooph <- c(unlist(stoc_run1_n2_28Aug[[10]]),unlist(stoc_run2_n10_28Aug[[10]]))

boxplot(stoc_R_burt, stoc_R_comp, stoc_D_div, stoc_A_del, stoc_A_fis, stoc_A_fra, stoc_C_spis, stoc_C_sta, stoc_Dicro, stoc_Ooph)

points(c(mean(stoc_R_burt), mean(stoc_R_comp), 
         mean(stoc_D_div), mean(stoc_A_del), 
         mean(stoc_A_fis), mean(stoc_A_fra), 
         mean(stoc_C_spis), mean(stoc_C_sta, na.rm = T),
         mean(stoc_Dicro), mean(stoc_Ooph)),
       pch = 19,
       col = "red")





#Average species determinism ####


example_mcmc_control <- list(n.burnin = 10, n.iteration = 1000, 
                             n.thin = 2)
auc_burt <- list()
auc_comp <- list()
auc_div <- list()
auc_del <- list()
auc_fis <- list()
auc_fram <- list()
auc_spis <- list()
auc_stam <- list()
auc_dic <- list()
auc_ooph <- list()
sp_num <- list()

cores <- detectCores()
cl<-makeCluster(cores - 1)
registerDoParallel(cl)

num <- 1
repeat {

  #Create datasubset  
  set.seed(Sys.time())
  rplots <- sample(150, 100)
  train <- all_dat[rplots, ]
  test <-  all_dat[-c(rplots), ]
  
  #Skip loop interation
  
  if(train %>% summarise(C_staminodiosum = sum(C_staminodiosum)) == 10) {
    next
  }
  
  #Record the number of species used in each plot
  sp_num[[num]] <- train  %>% summarise(R_burtoniae = sum(R_burtoniae),
                                        R_comptonii = sum(R_comptonii),
                                        D_diversifolium = sum(D_diversifolium),
                                        A_delaetii = sum(A_delaetii),
                                        A_fissum = sum(A_fissum),
                                        A_framesii = sum(A_framesii),
                                        C_spissum = sum(C_spissum),
                                        C_staminodiosum = sum(C_staminodiosum),
                                        Dicrocaulon_sp = sum(Dicrocaulon_sp),
                                        Oophytum_sp = sum(Oophytum_sp)
  ) 
  

  
  
  sp <- as.matrix(train[,1:10])
  covar <- as.matrix(train[,13:length(train)])
  
  #Run boral model
  mod <- boral(
    sp,
    covar,
    num.lv = 3,
    family = "binomial",
    #mcmc.control = example_mcmc_control,
    save.model = T
  )
  
  #Predict to test data
  test_covar <-  as.matrix(test[,13:length(test)])
  
  newpred <- predict.boral(mod, 
                           newX = test_covar, 
                           predict.type = "marginal",
                           est = "mean") 
  
  #Calculate auc for each species
  auc_burt[[num]] <- pROC::roc(test$R_burtoniae, newpred$linpred[,1]) %>% pROC::auc()
  auc_comp[[num]] <- pROC::roc(test$R_comptonii, newpred$linpred[,2]) %>% pROC::auc()
  auc_div[[num]] <- pROC::roc(test$D_diversifolium, newpred$linpred[,3]) %>% pROC::auc()
  auc_del[[num]] <- pROC::roc(test$A_delaetii, newpred$linpred[,4]) %>% pROC::auc()
  auc_fis[[num]] <- pROC::roc(test$A_fissum, newpred$linpred[,5]) %>% pROC::auc()
  auc_fram[[num]] <- pROC::roc(test$A_framesii, newpred$linpred[,6]) %>% pROC::auc()
  auc_spis[[num]] <- pROC::roc(test$C_spissum, newpred$linpred[,7]) %>% pROC::auc()
  auc_stam[[num]] <- pROC::roc(test$C_staminodiosum, newpred$linpred[,8]) %>% pROC::auc()
  auc_dic[[num]] <- pROC::roc(test$Dicrocaulon_sp, newpred$linpred[,9]) %>% pROC::auc()
  auc_ooph[[num]] <- pROC::roc(test$Oophytum_sp, newpred$linpred[,10]) %>% pROC::auc()
  
 
  
  num <- num + 1
  if(num > 40) {
    break
  }
}


#Stop parallel
stopCluster(cl)
registerDoSEQ()

##auc_scale32_randplots_28Aug <- list(auc_burt,
#                              auc_comp,
#                              auc_div ,
#                              auc_del ,
#                              auc_fis ,
#                              auc_fram,
#                              auc_spis,
#                              auc_stam,
#                              auc_dic ,
#                              auc_ooph)
#save(auc_scale32_randplots_28Aug, 
#     file = "/Users/larawootton/Documents/Honours/Data/auc_scale23_randplots_28Aug.rda")
#save(auc_scale10_run1_28Aug, 
         file = "/Users/larawootton/Documents/Honours/Data/auc_scale10_run1_28Aug.rda")
#save(auc_scale10_run2_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/auc_scale10_run2_28Aug.rda")
save(auc_scale24_run3_28Aug, file =  "/Users/larawootton/Documents/Honours/Data/auc_scale24_run3_28Aug.rda")
#save(sp_num_scale10_run1_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/sp_num_scale10_run1_28Aug.rda")
#save(sp_num_scale10_run2_28Aug, 
     file = "/Users/larawootton/Documents/Honours/Data/sp_num_scale10_run2_28Aug.rda")
save(sp_num_scale24_run3_28Aug, file = "/Users/larawootton/Documents/Honours/Data/sp_num_scale24_run3_28Aug.rda")
load("/Users/larawootton/Documents/Honours/Data/auc_scale10_run1_28Aug.rda")
load("/Users/larawootton/Documents/Honours/Data/auc_scale10_run2_28Aug.rda")
load("/Users/larawootton/Documents/Honours/Data/auc_scale24_run3_28Aug.rda")

load("/Users/larawootton/Documents/Honours/Data/sp_num_scale10_run1_28Aug.rda")
load("/Users/larawootton/Documents/Honours/Data/sp_num_scale10_run2_28Aug.rda")
load("/Users/larawootton/Documents/Honours/Data/sp_num_scale24_run3_28Aug.rda")

#auc_scale10_run1_28Aug<- list(auc_burt,
                              auc_comp,
                              auc_div ,
                              auc_del ,
                              auc_fis ,
                              auc_fram,
                              auc_spis,
                              auc_stam,
                              auc_dic ,
                              auc_ooph)
#auc_scale10_run2_28Aug<- list(auc_burt,
                              auc_comp,
                              auc_div ,
                              auc_del ,
                              auc_fis ,
                              auc_fram,
                              auc_spis,
                              auc_stam,
                              auc_dic ,
                              auc_ooph)
auc_scale24_run3_28Aug<- list(auc_burt,
                              auc_comp,
                              auc_div ,
                              auc_del ,
                              auc_fis ,
                              auc_fram,
                              auc_spis,
                              auc_stam,
                              auc_dic ,
                              auc_ooph)
#sp_num_scale10_run1_28Aug <- sp_num
#sp_num_scale10_run2_28Aug <- sp_num
sp_num_scale24_run3_28Aug <- sp_num


R_burt <- c(unlist(auc_scale10_run1_28Aug[[1]]), unlist(auc_scale10_run2_28Aug[[1]]), unlist(auc_scale24_run3_28Aug[[1]]))
R_comp <- c(unlist(auc_scale10_run1_28Aug[[2]]), unlist(auc_scale10_run2_28Aug[[2]]), unlist(auc_scale24_run3_28Aug[[2]]))
D_div <-  c(unlist(auc_scale10_run1_28Aug[[3]]), unlist(auc_scale10_run2_28Aug[[3]]), unlist(auc_scale24_run3_28Aug[[3]]))
A_del <-  c(unlist(auc_scale10_run1_28Aug[[4]]), unlist(auc_scale10_run2_28Aug[[4]]), unlist(auc_scale24_run3_28Aug[[4]]))
A_fis <-  c(unlist(auc_scale10_run1_28Aug[[5]]), unlist(auc_scale10_run2_28Aug[[5]]), unlist(auc_scale24_run3_28Aug[[5]]))
A_fra <-  c(unlist(auc_scale10_run1_28Aug[[6]]), unlist(auc_scale10_run2_28Aug[[6]]), unlist(auc_scale24_run3_28Aug[[6]]))
C_spis <- c(unlist(auc_scale10_run1_28Aug[[7]]), unlist(auc_scale10_run2_28Aug[[7]]), unlist(auc_scale24_run3_28Aug[[7]]))
C_sta <-  c(unlist(auc_scale10_run1_28Aug[[8]]), unlist(auc_scale10_run2_28Aug[[8]]), unlist(auc_scale24_run3_28Aug[[8]]))
Dicro <-  c(unlist(auc_scale10_run1_28Aug[[9]]), unlist(auc_scale10_run2_28Aug[[9]]), unlist(auc_scale24_run3_28Aug[[9]]))
Ooph <-   c(unlist(auc_scale10_run1_28Aug[[10]]),unlist(auc_scale10_run2_28Aug[[10]]), unlist(auc_scale24_run3_28Aug[[10]]))

sp_R_burt <- c(unlist(list.select(sp_num_scale10_run1_28Aug, R_burtoniae)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], R_burtoniae)),unlist(list.select(sp_num_scale24_run3_28Aug, R_burtoniae)))
sp_R_comp <- c(unlist(list.select(sp_num_scale10_run1_28Aug, R_comptonii)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], R_comptonii)),unlist(list.select(sp_num_scale24_run3_28Aug, R_comptonii)))
sp_D_div <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, D_diversifolium)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], D_diversifolium)),unlist(list.select(sp_num_scale24_run3_28Aug, D_diversifolium)))
sp_A_del <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, A_delaetii)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], A_delaetii)),unlist(list.select(sp_num_scale24_run3_28Aug, A_delaetii)))
sp_A_fis <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, A_fissum)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], A_fissum)),unlist(list.select(sp_num_scale24_run3_28Aug[1:23], A_fissum)))
sp_A_fra <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, A_framesii)),
              unlist(list.select(sp_num_scale10_run2_28Aug[1:10], A_framesii)),unlist(list.select(sp_num_scale24_run3_28Aug, A_framesii)))
sp_C_spis <- c(unlist(list.select(sp_num_scale10_run1_28Aug, C_spissum)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], C_spissum)),unlist(list.select(sp_num_scale24_run3_28Aug, C_spissum)))
sp_C_sta <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, C_staminodiosum)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], C_staminodiosum)),unlist(list.select(sp_num_scale24_run3_28Aug, C_staminodiosum)))
sp_Dicro <-  c(unlist(list.select(sp_num_scale10_run1_28Aug, Dicrocaulon_sp)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], Dicrocaulon_sp)),unlist(list.select(sp_num_scale24_run3_28Aug, Dicrocaulon_sp)))
sp_Ooph <-   c(unlist(list.select(sp_num_scale10_run1_28Aug, Oophytum_sp)),
               unlist(list.select(sp_num_scale10_run2_28Aug[1:10], Oophytum_sp)),unlist(list.select(sp_num_scale24_run3_28Aug,  Oophytum_sp)))

plot(R_burt ~ sp_R_burt,
     xlab = "R. burtoniae sample size",
     ylab = "AUC",
     pch = 19)
modb <- lm(R_burt ~ sp_R_burt)
summary(modb)
abline(modb)

plot(R_comp ~ sp_R_comp,
     xlab = "R. comptonii sample size",
     ylab = "AUC",
     pch = 19)

modrc <- lm(R_comp ~ sp_R_comp)
summary(modrc)


plot(D_div ~ sp_D_div,
     xlab = "D. diversifolium sample size",
     ylab = "AUC",
     pch = 19)
modd <- lm(D_div ~ sp_D_div)
summary(modd)

plot(A_del ~ sp_A_del,
     xlab = "A. delaetii sample size",
     ylab = "AUC",
     pch = 19)
modad <- lm(A_del ~ sp_A_del)
summary(modad)


plot(A_fis ~ sp_A_fis,
     xlab = "A. fissum sample size",
     ylab = "AUC",
     pch = 19)
modaf <- lm(A_fis ~ sp_A_fis)
summary(modaf)

plot(A_fra ~ sp_A_fra,
     xlab = "A. framesii sample size",
     ylab = "AUC",
     pch = 19)
moda <- lm(A_fra ~ sp_A_fra)
summary(moda)
abline(moda)


plot(C_spis ~ sp_C_spis,
     xlab = "C. spissum sample size",
     ylab = "AUC",
     pch = 19)
modc <- lm(C_spis ~ sp_C_spis)
summary(modc)

plot(C_sta ~ sp_C_sta,
     xlab = "C. staminodiosum sample size",
     ylab = "AUC",
     pch = 19)
modcs <- lm(C_sta ~ sp_C_sta)
summary(modcs)

plot(Dicro ~ sp_Dicro,
     xlab = "Dicrocaulon sample size",
     ylab = "AUC",
     pch = 19)
moddi <- lm(Dicro ~ sp_Dicro)
summary(moddi)

plot(Ooph ~ sp_Ooph,
     xlab = "Oophytum sample size",
     ylab = "AUC",
     pch = 19)
modo <- lm(Ooph ~ sp_Ooph)
summary(modo)



boxplot(R_burt, R_comp, D_div, A_del, A_fis, A_fra, C_spis, C_sta, Dicro, Ooph)
points(c(mean(R_burt), mean(R_comp), 
       mean(D_div), mean(A_del), 
       mean(A_fis), mean(A_fra), 
       mean(C_spis), mean(C_sta),
       mean(Dicro), mean(Ooph)),
       pch = 19,
       col = "red")

boxplot(R_burt, stoc_R_burt, 
        R_comp, stoc_R_comp,
        D_div, stoc_D_div,
        A_del, stoc_A_del,
        A_fis, stoc_A_fis,
        A_fra, stoc_A_fra,
        C_spis, stoc_C_spis,
        C_sta, stoc_C_sta,
        Dicro, stoc_Dicro,
        Ooph, stoc_Ooph,
        col = c("white", "grey"))

#Outlying means index (OMI) ####

# Measures the distance between the mean habitat conditions used by each species 
#and the mean habitat conditions of the study area


#First make a pca of the soil variables

vars <- all_dat %>% select(13:26)

pca_soil <- dudi.pca(vars, scannf = FALSE)
scatter(pca_soil)

#Percentage variation associated with each axis
100 * pca_soil$eig/sum(pca_soil$eig)

#Correlation between each variable and the axes
pca_soil$co

nic <- canomi(pca_soil, all_dat[,1:10], scannf=FALSE)
nic
plot(nic)


omi <- niche(pca_soil, all_dat[,1:10], scannf = FALSE)
plot(omi)

#variation of the relationship between species and environmental gradients
100 * omi$eig/sum(omi$eig)


par(mfrow=c(1,1))
sco.distri(omi$ls[,1],all_dat[,1:10],clab=0.7)
sco.distri(omi$ls[,2],all_dat[,1:10],clab=0.7)


s.corcircle(omi$as, 1, 2, sub = "Axis", csub = 2, 
            clabel = 1.25)
s.arrow(omi$c1, 1, 2, sub = "Variables", csub = 2, 
        clabel = 1.25)
scatterutil.eigen(omi$eig, wsel = c(1, 2))
s.label(omi$ls, 1, 2, clabel = 0, cpoint = 2, sub = "Samples and Species", 
        csub = 2)
s.label(omi$li, 1, 2, clabel = 1.5, add.plot = TRUE)
s.label(omi$ls, 1, 2, clabel = 1.25, sub = "Samples", 
        csub = 2)
s.distri(omi$ls, eval.parent(as.list(omi$call)[[3]]), 
         cstar = 0, axesell = FALSE, cellipse = 1, sub = "Niches", csub = 2)
s.label(omi$li, 1, 2, clabel = 0.7, add.plot = TRUE)
