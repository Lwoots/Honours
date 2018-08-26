#Interrogating prediction differences
#Cause of discrepancy between predictive ability
#Honours 2018

#Set up -----------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(ggplot2, MASS, dplyr, boral, pROC)

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

dat <- left_join(sp_occ, my_soil_df[,c(1, 5)], by = "plot") %>% 
  filter(!is.na(type),!plot == "R40 ")

all_dat <- left_join(sp_occ, my_soil_df %>% select(colnames(soil_dat[2:15]), plot), by = "plot") %>% 
  na.omit()



#Do the sites have different soil properties?
#Linear discriminant analysis
#LDA####

lda_mod <- lda(site ~ ., data = soil_dat)
cols <- c("black", "red", "green")[soil_dat$site]
plot(lda_mod, col = cols)


#How does the number of occurrences per species differ between sites? ####

total <- apply(X = dat[,1:10], MARGIN = 2, FUN = sum)

site_sums <- dat %>% group_by(site) %>% summarise(R_burtoniae = sum(R_burtoniae),
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
site_sums <- rbind(site_sums[,2:11], tot = total)

sum_dat <- as.data.frame(t(site_sums[,1:10]))
colnames(sum_dat) <- c("Site 1", "Site 2", "Site 3", "Total")
sum_dat


#How deterministic are the species? ####

#Practice loop

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



num <- 1
repeat {
  set.seed(Sys.time())
  rplots <- sample(150, 100)
  train <- all_dat[rplots, ]
  test <-  all_dat[-c(rplots), ]
  
  sp <- as.matrix(train[,1:10])
  covar <- as.matrix(train[,13:length(train)])
  
  mod <- boral(
    sp,
    covar,
    num.lv = 3,
    family = "binomial",
    #mcmc.control = example_mcmc_control,
    save.model = T
  )
  
  test_covar <-  as.matrix(test[,13:length(test)])
  
  newpred <- predict.boral(mod, 
                           newX = test_covar, 
                           predict.type = "marginal",
                           est = "mean") 
 
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
  if(num > 4) {
    break
  }
}


bl <- train %>% summarise(R_burtoniae = sum(R_burtoniae),
          R_comptonii = sum(R_comptonii),
          D_diversifolium = sum(D_diversifolium),
          A_delaetii = sum(A_delaetii),
          A_fissum = sum(A_fissum),
          A_framesii = sum(A_framesii),
          C_spissum = sum(C_spissum),
          C_staminodiosum = sum(C_staminodiosum),
          Dicrocaulon_sp = sum(Dicrocaulon_sp),
          Oophytum_sp = sum(Oophytum_sp))
bt <- test %>% summarise(R_burtoniae = sum(R_burtoniae),
                                   R_comptonii = sum(R_comptonii),
                                   D_diversifolium = sum(D_diversifolium),
                                   A_delaetii = sum(A_delaetii),
                                   A_fissum = sum(A_fissum),
                                   A_framesii = sum(A_framesii),
                                   C_spissum = sum(C_spissum),
                                   C_staminodiosum = sum(C_staminodiosum),
                                   Dicrocaulon_sp = sum(Dicrocaulon_sp),
                                   Oophytum_sp = sum(Oophytum_sp))

hist(unlist(auc_comp))


newpred$linpred[,1] %>% pROC::roc(test$R_comptonii, newpred$linpred[,2]) %>% pROC::auc()

pROC::roc(test$R_comptonii, newpred$linpred[,2]) %>% pROC::auc()

comp <- prediction(newpred$linpred[,2], test$R_comptonii) %>% performance("auc")




aucs@y.values
ls <- list(unlist(aucs@y.values), unlist(aucs@y.values), c(8,9,40))

ls[[2]]
ls[1:3][3]
t <- unlist(lapply(ls, "[[", 1))
hist(t)

rank(t)

ls[[1]][[1]]


roc()

