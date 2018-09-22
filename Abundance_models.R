#More abundance models

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

all_dat <- left_join(nb_sp, my_soil_df %>% select(colnames(soil_dat[2:16]), plot), by = "plot") %>% 
  na.omit()
all_dat[,13:27] <- scale(all_dat[,13:27])

#Standardise abundance by dividing by the sd of plots not including zero
#all_dat$R_burtoniae <- all_dat$R_burtoniae/(all_dat$R_burtoniae[which(all_dat$R_burtoniae > 0)] %>% sd())
#all_dat$R_comptonii <- all_dat$R_comptonii/(all_dat$R_comptonii[which(all_dat$R_comptonii > 0)] %>% sd())
#all_dat$D_diversifolium <- all_dat$D_diversifolium/(all_dat$D_diversifolium[which(all_dat$D_diversifolium > 0)] %>% #sd())
#all_dat$A_delaetii <- all_dat$A_delaetii/(all_dat$A_delaetii[which(all_dat$A_delaetii > 0)] %>% sd())
#all_dat$A_fissum <- all_dat$A_fissum/(all_dat$A_fissum[which(all_dat$A_fissum > 0)] %>% sd())
#all_dat$A_framesii <- all_dat$A_framesii/(all_dat$A_framesii[which(all_dat$A_framesii > 0)] %>% sd())
#all_dat$C_spissum <- all_dat$C_spissum/(all_dat$C_spissum[which(all_dat$C_spissum > 0)] %>% sd())
#all_dat$C_staminodiosum <- all_dat$C_staminodiosum/(all_dat$C_staminodiosum[which(all_dat$C_staminodiosum > 0)] %>% #sd())
#all_dat$Dicrocaulon_sp <- all_dat$Dicrocaulon_sp/(all_dat$Dicrocaulon_sp[which(all_dat$Dicrocaulon_sp > 0)] %>% sd())
#all_dat$Oophytum_sp <- all_dat$Oophytum_sp/(all_dat$Oophytum_sp[which(all_dat$Oophytum_sp > 0)] %>% sd())

#Abundance model - normal ####

sp <- as.matrix(all_dat[,3:12])
covar <- as.matrix(all_dat[,13:length(all_dat)])

normal_adn_1Sept <- boral(sp,
                         X = covar,
                         family = "lnormal",
                         num.lv = 3,
                         save.model = T,
                         calc.ics = T)
save(normal_adn_1Sept, file = "normal_adn_1Sept.rda")
plot.boral(normal_adn_1Sept)

norm_fit <- fitted.boral(normal_adn_1Sept, est = "mean")

pred <- as.data.frame(norm_fit$out)

plot(pred$R_burtoniae ~ all_dat$R_burtoniae,
     ylab = "Predicted",
     xlab = "",
     main = "R_burtoniae")
abline(0,1)
plot(pred$R_comptonii ~ all_dat$R_comptonii,
     xlab = "",
     ylab = "",
     main = "R_comptonii")
abline(0,1)
plot(pred$D_diversifolium ~ all_dat$D_diversifolium,
     ylab = "Predicted",
     xlab = "",
     main = "D_diversifolium")
abline(0,1)
plot(pred$A_delaetii ~ all_dat$A_delaetii,
     xlab = "",
     ylab = "",
     main = "A_delaetii")
abline(0,1)
plot(pred$A_fissum ~ all_dat$A_fissum,
     ylab = "Predicted",
     xlab = "",
     main = "A_fissum")
abline(0,1)
plot(pred$A_framesii ~ all_dat$A_framesii,
     xlab = "",
     ylab = "",
     main = "A_framesii")
abline(0,1)
plot(pred$C_spissum ~ all_dat$C_spissum,
     ylab = "Predicted",
     xlab = "",
     main = "C_spissum")
abline(0,1)
plot(pred$C_staminodiosum ~ all_dat$C_staminodiosum,
     xlab = "",
     ylab = "",
     main = "C_staminodiosum")
abline(0,1)
plot(pred$Dicrocaulon_sp ~ all_dat$Dicrocaulon_sp,
     ylab = "Predicted",
     xlab = "Observed",
     main = "Dicrocaulon_sp")
abline(0,1)
plot(pred$Oophytum_sp ~ all_dat$Oophytum_sp,
     ylab = "Predicted",
     xlab = "",
     main = "Oophytum_sp")
abline(0,1)


#Negative binomial ####

sp <- as.matrix(all_dat[,3:12])
covar <- as.matrix(all_dat[,13:length(all_dat)])

negbio_abn_1Sept <- boral(sp,
                          X = covar,
                          family = "negative.binomial",
                          num.lv = 3,
                          save.model = T,
                          calc.ics = T)
save(negbio_abn_1Sept, file = "negbio_abn_1Sept.rda")
load("negbio_abn_1Sept.rda")
plot.boral(negbio_abn_1Sept)
abline(0,1)

var <- calc.varpart(negbio_abn_1Sept)
norm_fit <- fitted.boral(negbio_abn_1Sept, est = "mean")

pred <- as.data.frame(norm_fit$out) 

plot(pred$R_burtoniae ~ all_dat$R_burtoniae,
     ylab = "Predicted",
     xlab = "",
     main = "R_burtoniae")
abline(0,1)
plot(pred$R_comptonii ~ all_dat$R_comptonii,
     xlab = "",
     ylab = "",
     main = "R_comptonii")
abline(0,1)
plot(pred$D_diversifolium ~ all_dat$D_diversifolium,
     ylab = "Predicted",
     xlab = "",
     main = "D_diversifolium")
abline(0,1)
plot(pred$A_delaetii ~ all_dat$A_delaetii,
     xlab = "",
     ylab = "",
     main = "A_delaetii")
abline(0,1)
plot(pred$A_fissum ~ all_dat$A_fissum,
     ylab = "Predicted",
     xlab = "",
     main = "A_fissum")
abline(0,1)
plot(pred$A_framesii ~ all_dat$A_framesii,
     xlab = "",
     ylab = "",
     main = "A_framesii")
abline(0,1)
plot(pred$C_spissum ~ all_dat$C_spissum,
     ylab = "Predicted",
     xlab = "",
     main = "C_spissum")
abline(0,1)
plot(pred$C_staminodiosum ~ all_dat$C_staminodiosum,
     xlab = "",
     ylab = "",
     main = "C_staminodiosum")
abline(0,1)
plot(pred$Dicrocaulon_sp ~ all_dat$Dicrocaulon_sp,
     ylab = "Predicted",
     xlab = "Observed",
     main = "Dicrocaulon_sp")
abline(0,1)
plot(pred$Oophytum_sp ~ all_dat$Oophytum_sp,
     ylab = "Predicted",
     xlab = "",
     main = "Oophytum_sp")
abline(0,1)


#Poisson ####

sp <- as.matrix(all_dat[,3:12])
covar <- as.matrix(all_dat[,13:length(all_dat)])

poisson_abn_1Sept <- boral(sp,
                          X = covar,
                          family = "poisson",
                          num.lv = 3,
                          save.model = T,
                          calc.ics = T)
save(poisson_abn_1Sept, file = "poisson_abn_1Sept.rda")
plot.boral(poisson_abn_1Sept)

pois_fit <- fitted.boral(poisson_abn_1Sept, est = "mean")

pred <- as.data.frame(pois_fit$out) 

plot(pred$R_burtoniae ~ all_dat$R_burtoniae,
     ylab = "Predicted",
     xlab = "",
     main = "R_burtoniae")
abline(0,1)
plot(pred$R_comptonii ~ all_dat$R_comptonii,
     xlab = "",
     ylab = "",
     main = "R_comptonii")
abline(0,1)
plot(pred$D_diversifolium ~ all_dat$D_diversifolium,
     ylab = "Predicted",
     xlab = "",
     main = "D_diversifolium")
abline(0,1)
plot(pred$A_delaetii ~ all_dat$A_delaetii,
     xlab = "",
     ylab = "",
     main = "A_delaetii")
abline(0,1)
plot(pred$A_fissum ~ all_dat$A_fissum,
     ylab = "Predicted",
     xlab = "",
     main = "A_fissum")
abline(0,1)
plot(pred$A_framesii ~ all_dat$A_framesii,
     xlab = "",
     ylab = "",
     main = "A_framesii")
abline(0,1)
plot(pred$C_spissum ~ all_dat$C_spissum,
     ylab = "Predicted",
     xlab = "",
     main = "C_spissum")
abline(0,1)
plot(pred$C_staminodiosum ~ all_dat$C_staminodiosum,
     xlab = "",
     ylab = "",
     main = "C_staminodiosum")
abline(0,1)
plot(pred$Dicrocaulon_sp ~ all_dat$Dicrocaulon_sp,
     ylab = "Predicted",
     xlab = "Observed",
     main = "Dicrocaulon_sp")
abline(0,1)
plot(pred$Oophytum_sp ~ all_dat$Oophytum_sp,
     ylab = "Predicted",
     xlab = "",
     main = "Oophytum_sp")
abline(0,1)

#tweedie ####

sp <- as.matrix(all_dat[,3:12])
covar <- as.matrix(all_dat[,13:length(all_dat)])

twee_20sept <- boral(sp,
                          X = covar,
                          family = "exponential",
                          num.lv = 3,
                          save.model = T,
                          calc.ics = T)

plot.boral(twee_20sept)
abline(0,1)


#Randomisations ####


rplots <- sample(150, 100)
train <- all_dat[rplots, 3:27]
test <-  all_dat[-c(rplots),3:27]

train <- all_dat %>% filter(site == "site3")
test <- all_dat %>% filter(site == "site2")

sp <- train[,3:12]
covar <- train[,13:27]



mod <- boral(sp,
                   X = covar,
                   family = "poisson",
                   num.lv = 3,
                   save.model = T)

test_covar <-  as.matrix(test[,13:length(test)])

newpred <- predict.boral(mod, 
                         newX = test_covar, 
                         predict.type = "marginal",
                         est = "mean")
newpred$linpred <- exp(newpred$linpred)

plot(newpred$linpred[,8] ~ test$C_staminodiosum)

t <- newpred$linpred[,1] - test[,3]
t <- ( test[,6] - newpred$linpred[,3] ) %>% abs() %>% mean()

mean(t)
mean(abs(t))

sum(test[,1])

#The reps ####

diff_burt <- list()
diff_comp <- list()
diff_div <- list()
diff_del <- list()
diff_fis <- list()
diff_fra <- list()
diff_spi <- list()
diff_sta <- list()
diff_dic <- list()
diff_oop <- list()
sp_num <- list()

num <- 1

repeat {
  
  #Create data subset  
  set.seed(Sys.time())
  rplots <- sample(150, 100)
  train <- all_dat[rplots, 3:27]
  test <-  all_dat[-c(rplots), 3:27]
 
  #Record the number of individuals per species used in each model
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
  
  #Set up boral
  
  sp <- train[,1:10]
  covar <- train[,11:length(train)]
  
  mod <- boral(sp,
               X = covar,
               family = "poisson",
               num.lv = 3,
               save.model = T)
  
  #Make predictions
  test_covar <- test[,11:length(test)]
  
  newpred <- predict.boral(mod, 
                           newX = test_covar, 
                           predict.type = "marginal",
                           est = "mean")
  
  newpred$linpred <- exp(newpred$linpred)
  
  #Record mean difference between 
  
  diff_burt[[num]] <- (newpred$linpred[,1] - test[,1]) %>% abs() %>% mean()
  diff_comp[[num]] <- (newpred$linpred[,2] - test[,2]) %>% abs() %>% mean()
  diff_div[[num]]  <- (newpred$linpred[,3] - test[,3]) %>% abs() %>% mean()
  diff_del[[num]]  <- (newpred$linpred[,4] - test[,4]) %>% abs() %>% mean()
  diff_fis[[num]]  <- (newpred$linpred[,5] - test[,5]) %>% abs() %>% mean()
  diff_fra[[num]]  <- (newpred$linpred[,6] - test[,6]) %>% abs() %>% mean()
  diff_spi[[num]]  <- (newpred$linpred[,7] - test[,7]) %>% abs() %>% mean()
  diff_sta[[num]]  <- (newpred$linpred[,8] - test[,8]) %>% abs() %>% mean()
  diff_dic[[num]]  <- (newpred$linpred[,9] - test[,9]) %>% abs() %>% mean()
  diff_oop[[num]]  <- (newpred$linpred[,10] - test[,10]) %>% abs() %>% mean()
  
  num <- num + 1
  if(num > 54) {
    break
  }
}

diff_run2_n54_7Sept <- list(diff_burt,
                            diff_comp,
                            diff_div, 
                            diff_del, 
                            diff_fis, 
                            diff_fra, 
                            diff_spi, 
                            diff_sta, 
                            diff_dic, 
                            diff_oop )

save(diff_run1_n45_6Sept, file = "diff_run1_n45_6Sept.rda")
save(diff_run2_n54_7Sept, file = "diff_run2_n54_7Sept.rda")
load("diff_run1_n45_6Sept.rda")
load("diff_run2_n54_7Sept.rda")

sp_diff_run1_n45_6Sept <- sp_num
sp_diff_run2_n54_7Sept <- sp_num
save(sp_diff_run1_n45_6Sept, file = "sp_diff_run1_n45_6Sept.rda")
save(sp_diff_run2_n54_7Sept, file = "sp_diff_run2_n54_7Sept.rda")
load("sp_diff_run1_n45_6Sept.rda")
load("sp_diff_run2_n54_7Sept.rda")

#Extract differences
R_bur <- c(unlist(diff_run1_n45_6Sept[[1]]), unlist(diff_run2_n54_7Sept[[1]]))
R_com <- c(unlist(diff_run1_n45_6Sept[[2]]), unlist(diff_run2_n54_7Sept[[2]])) 
D_div <- c(unlist(diff_run1_n45_6Sept[[3]]), unlist(diff_run2_n54_7Sept[[3]]))
A_del <- c(unlist(diff_run1_n45_6Sept[[4]]), unlist(diff_run2_n54_7Sept[[4]]))
A_fis <- c(unlist(diff_run1_n45_6Sept[[5]]), unlist(diff_run2_n54_7Sept[[5]]))
A_fra <- c(unlist(diff_run1_n45_6Sept[[6]]), unlist(diff_run2_n54_7Sept[[6]]))
C_spi <- c(unlist(diff_run1_n45_6Sept[[7]]), unlist(diff_run2_n54_7Sept[[7]]))
C_sta <- c(unlist(diff_run1_n45_6Sept[[8]]), unlist(diff_run2_n54_7Sept[[8]]))
Dicro <- c(unlist(diff_run1_n45_6Sept[[9]]), unlist(diff_run2_n54_7Sept[[9]]))
Ooph <- c(unlist(diff_run1_n45_6Sept[[10]]), unlist(diff_run2_n54_7Sept[[10]]))

sp_R_bur <- c(unlist(list.select(sp_diff_run1_n45_6Sept, R_burtoniae)), unlist(list.select(sp_diff_run2_n54_7Sept, R_burtoniae)))
sp_R_com <- c(unlist(list.select(sp_diff_run1_n45_6Sept, R_comptonii)))
sp_D_div <- c(unlist(list.select(sp_diff_run1_n45_6Sept, D_diversifolium)))
sp_A_del <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_delaetii)))
sp_A_fis <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_fissum)))
sp_A_fra <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_framesii)))
sp_C_spi <- c(unlist(list.select(sp_diff_run1_n45_6Sept, C_spissum)))
sp_C_sta <- c(unlist(list.select(sp_diff_run1_n45_6Sept, C_staminodiosum)))
sp_Dicro <- c(unlist(list.select(sp_diff_run1_n45_6Sept, Dicrocaulon_sp)))
sp_Ooph <- c(unlist(list.select(sp_diff_run1_n45_6Sept, Oophytum_sp)))

plot(sp_R_bur, R_bur)



boxplot(R_bur, 
        R_com, 
        D_div, 
        #A_del, 
        #A_fis, 
        A_fra, 
        C_spi, 
        C_sta, 
        Dicro, 
        Ooph)

boxplot(log(A_del), log(A_fra),  
        log(C_spi), log(Dicro), 
        log(Ooph), log(D_div), 
        log(R_com), log(R_bur), 
        log(A_fis), log(C_sta))
abline(1, 0)


mean(log(R_bur))
mean(log(R_com))
mean(log(D_div))
mean(log(A_del))
mean(log(A_fis))
mean(log(A_fra))
mean(log(C_spi))
mean(log(C_sta))
mean(log(Dicro))
mean(log(Ooph))



#Stochastic null model ####

stoc_diff_bur <- list()
stoc_diff_com <- list()
stoc_diff_div <- list()
stoc_diff_del <- list()
stoc_diff_fis <- list()
stoc_diff_fra <- list()
stoc_diff_spi <- list()
stoc_diff_sta <- list()
stoc_diff_dic <- list()
stoc_diff_oop <- list()
stoc_sp_num <- list()

num = 1

repeat {
  
  #Shuffle species, keeping row and column sums equal to original data.
  rsp <- permatswap(all_dat[,3:12], 
                    times = 1, 
                    burnin = 5, 
                    fixedmar = "both", 
                    mtype = "count")
  rand_dat <- cbind(rsp$perm, all_dat[,13:27])
  
  #Choose subset of plots
  set.seed(Sys.time())
  rplots <- sample(150, 100)
  train <- rand_dat[rplots,]
  test <- rand_dat[-c(rplots),]
  
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
  
  #Set up JSDM
  sp <- as.matrix(train[,1:10])
  covar <- as.matrix(train[,11:length(train)])
  
  mod <- boral(sp,
               X = covar,
               family = "poisson",
               num.lv = 3,
               save.model = T)
  
  #Make predictions
  test_covar <- test[,11:length(test)]
  
  newpred <- predict.boral(mod, 
                           newX = test_covar, 
                           predict.type = "marginal",
                           est = "mean")
  
  newpred$linpred <- exp(newpred$linpred)
  
  #Record mean difference between 
  
  stoc_diff_bur[[num]] <- (newpred$linpred[,1] - test[,1]) %>% abs() %>% mean()
  stoc_diff_com[[num]] <- (newpred$linpred[,2] - test[,2]) %>% abs() %>% mean()
  stoc_diff_div[[num]]  <- (newpred$linpred[,3] - test[,3]) %>% abs() %>% mean()
  stoc_diff_del[[num]]  <- (newpred$linpred[,4] - test[,4]) %>% abs() %>% mean()
  stoc_diff_fis[[num]]  <- (newpred$linpred[,5] - test[,5]) %>% abs() %>% mean()
  stoc_diff_fra[[num]]  <- (newpred$linpred[,6] - test[,6]) %>% abs() %>% mean()
  stoc_diff_spi[[num]]  <- (newpred$linpred[,7] - test[,7]) %>% abs() %>% mean()
  stoc_diff_sta[[num]]  <- (newpred$linpred[,8] - test[,8]) %>% abs() %>% mean()
  stoc_diff_dic[[num]]  <- (newpred$linpred[,9] - test[,9]) %>% abs() %>% mean()
  stoc_diff_oop[[num]]  <- (newpred$linpred[,10] - test[,10]) %>% abs() %>% mean()
  
  num <- num + 1
  if(num > 99) {
    break
  }
}


stoc_run1_n99_8Sept <- list(
  stoc_diff_bur,
  stoc_diff_com,
  stoc_diff_div,
  stoc_diff_del,
  stoc_diff_fis,
  stoc_diff_fra,
  stoc_diff_spi,
  stoc_diff_sta,
  stoc_diff_dic,
  stoc_diff_oop
)

save(stoc_run1_n99_8Sept, file = "stoc_run1_n99_8Sept.rda")
load("stoc_run1_n99_8Sept.rda")

stoc_R_bur <- unlist(stoc_run1_n99_8Sept[[1]])
stoc_R_com <- unlist(stoc_run1_n99_8Sept[[2]])
stoc_D_div <- unlist(stoc_run1_n99_8Sept[[3]])
stoc_A_del <- unlist(stoc_run1_n99_8Sept[[4]])
stoc_A_fis <- unlist(stoc_run1_n99_8Sept[[5]])
stoc_A_fra <- unlist(stoc_run1_n99_8Sept[[6]])
stoc_C_spi <- unlist(stoc_run1_n99_8Sept[[7]])
stoc_C_sta <- unlist(stoc_run1_n99_8Sept[[8]])
stoc_Dicro <- unlist(stoc_run1_n99_8Sept[[8]])
stoc_Ooph <- unlist(stoc_run1_n99_8Sept[[10]])



boxplot(log(stoc_R_bur), log(stoc_R_com),log(stoc_D_div),
        log(stoc_A_del), log(stoc_A_fis), log(stoc_A_fra), log(stoc_C_spi),
        log(stoc_C_sta), log(stoc_Dicro), log(stoc_Ooph))

median(stoc_R_bur)
mean(stoc_R_bur)

#Figure


boxplot(
        log(C_sta), log(stoc_C_sta),
        log(A_fis), log(stoc_A_fis),
        log(R_bur), log(stoc_R_bur),
        log(R_com), log(stoc_R_com),
        log(D_div), log(stoc_D_div),
        log(Ooph), log(stoc_Ooph),
        
        log(C_spi), log(stoc_C_spi),
        log(Dicro), log(stoc_Dicro),
        log(A_fra), log(stoc_A_fra),
        log(A_del), log(stoc_A_del),
        
        col = c("steelblue1", "grey"),
        whisklty = 1,
        staplelty = 0,
        at = c(1,2,4,5,7,8,10,11,13,14,16,17,19,20, 22,23,25,26,28,29),
        xaxt = "n",
        cex = 0.75,
        ylab = "Log( |predicted - observed counts| )")

axis(1, at = seq(1.5,28.5, 3),
     labels = F,
     tck = -0.01
)

labels <- c(
  "C. staminodiosum",
  "A. fissum",
  "R. burtoniae",
  "R. comptonii",
  "D. diversifolium",
  "Oophytum sp.",
  "C. spissum",
  "Dicrocaulon sp.",
  "A. framesii",
  "A. delaetii"
)

text(x = seq(1.5,28.5, 3), y = -5.8, srt = 40, adj= 1, xpd = TRUE, labels = labels, cex=0.9)
points(x = seq(1,28,3),
       y = sort(c(mean(log(R_bur)),
                  mean(log(R_com)), 
                  mean(log(D_div)), 
                  mean(log(A_del)), 
                  mean(log(A_fis)),
                  mean(log(A_fra)), 
                  mean(log(C_spi)), 
                  mean(log(C_sta)),
                  mean(log(Dicro)), 
                  mean(log(Ooph)))),
       pch = 21,
       col = "black",
       bg = "#F6B540",
       cex = 1.4
)

mtext(at = c(4.5, 7.5, 16.5, 19.5), 
      side = 3, 
      text = "*", 
      cex = 2)


t.test(log(C_sta), log(stoc_C_sta), alternative = "less")
t.test(log(A_fis), log(stoc_A_fis), alternative = "less") #*
t.test(log(R_bur), log(stoc_R_bur), alternative = "less") #*
t.test(log(R_com), log(stoc_R_com), alternative = "less")
t.test(log(D_div), log(stoc_D_div), alternative = "less")
t.test(log(Ooph), log(stoc_Ooph), alternative = "less")   #*
t.test(log(C_spi), log(stoc_C_spi), alternative = "less") #*
t.test(log(Dicro), log(stoc_Dicro), alternative = "less")
t.test(log(A_fra), log(stoc_A_fra), alternative = "less")
t.test(log(A_del), log(stoc_A_del), alternative = "less")



all_dat %>% filter(all_dat$C_staminodiosum > 0) %>% summarise(c = mean(C_staminodiosum))
all_dat %>% filter(all_dat$A_fissum > 0) %>% summarise(c = mean(A_fissum))
all_dat %>% filter(all_dat$R_burtoniae > 0) %>% summarise(c = mean(R_burtoniae))
all_dat %>% filter(all_dat$R_comptonii > 0) %>% summarise(c = mean(R_comptonii))
all_dat %>% filter(all_dat$D_diversifolium > 0) %>% summarise(c = mean(D_diversifolium))
all_dat %>% filter(all_dat$Oophytum_sp > 0) %>% summarise(c = mean(Oophytum_sp))
all_dat %>% filter(all_dat$C_spissum > 0) %>% summarise(c = mean(C_spissum))
all_dat %>% filter(all_dat$Dicrocaulon_sp > 0) %>% summarise(c = mean(Dicrocaulon_sp))
all_dat %>% filter(all_dat$A_framesii > 0) %>% summarise(c = mean(A_framesii))
all_dat %>% filter(all_dat$A_delaetii > 0) %>% summarise(c = mean(A_delaetii))


mean(all_dat$A_fissum)
mean(all_dat$R_burtoniae)
mean(all_dat$R_comptonii)
mean(all_dat$D_diversifolium)
mean(all_dat$Oophytum_sp)
mean(all_dat$C_spissum)
mean(all_dat$Dicrocaulon_sp)

plot(newpred$linpred[,5] ~ test$A_fissum)




C_sta_sd <- C_sta / (as.numeric(all_dat %>% filter(all_dat$C_staminodiosum > 0) %>% summarise(c = mean(C_staminodiosum))))
A_fis_sd <- A_fis / as.numeric(all_dat %>% filter(all_dat$A_fissum > 0) %>% summarise(c = mean(A_fissum)))
  


log(mean(A_fis))
max(A_fis)

mean(Ooph)
log(mean(R_com))
log(mean(A_del))


p <- permatswap(all_dat[,3:12], method = "swsh", fixedmar="columns", shuffle = "samp")

p$perm

c <- make.commsim()



#create a vector of random values
foo <- rnorm(n=100, mean=20, sd=5)
#randomly choose 15 indices to replace
#this is the step in which I thought I was clever
#because I use which() and %in% in the same line
ind <- which(foo %in% sample(foo, 15))
#now replace those indices in foo with NA
foo[ind]<-NA
#here is our vector with 15 random NAs 
foo


abn <- which(all_dat$R_burtoniae == 0 )

fo <- sample(all_dat$R_burtoniae[-abn], length(all_dat$R_burtoniae[-abn]))
shuf_burt <- vector(length = 150)
shuf_burt[which(all_dat$R_burtoniae == 0 )] <- 0
shuf_burt[which(all_dat$R_burtoniae > 0 )] <- fo
shuf_burt
