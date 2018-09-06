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



#Abundance model - normal ####

sp <- as.matrix(all_dat[,3:12])
covar <- as.matrix(all_dat[,13:length(all_dat)])

normal_adn_1Sept <- boral(sp,
                         X = covar,
                         family = "normal",
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
plot.boral(negbio_abn_1Sept)

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
  if(num > 45) {
    break
  }
}

diff_run1_n45_6Sept <- list(diff_burt,
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

sp_diff_run1_n45_6Sept <- sp_num
save(sp_diff_run1_n45_6Sept, file = "sp_diff_run1_n45_6Sept.rda")

#Extract differences
R_bur <- c(unlist(diff_run1_n45_6Sept[[1]]))
R_com <- c(unlist(diff_run1_n45_6Sept[[2]])) 
D_div <- c(unlist(diff_run1_n45_6Sept[[3]]))
A_del <- c(unlist(diff_run1_n45_6Sept[[4]]))
A_fis <- c(unlist(diff_run1_n45_6Sept[[5]]))
A_fra <- c(unlist(diff_run1_n45_6Sept[[6]]))
C_spi <- c(unlist(diff_run1_n45_6Sept[[7]]))
C_sta <- c(unlist(diff_run1_n45_6Sept[[8]]))
Dicro <- c(unlist(diff_run1_n45_6Sept[[9]]))
Ooph <- c(unlist(diff_run1_n45_6Sept[[10]]))

sp_R_bur <- c(unlist(list.select(sp_diff_run1_n45_6Sept, R_burtoniae)))
sp_R_com <- c(unlist(list.select(sp_diff_run1_n45_6Sept, R_comptonii)))
sp_D_div <- c(unlist(list.select(sp_diff_run1_n45_6Sept, D_diversifolium)))
sp_A_del <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_delaetii)))
sp_A_fis <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_fissum)))
sp_A_fra <- c(unlist(list.select(sp_diff_run1_n45_6Sept, A_framesii)))
sp_C_spi <- c(unlist(list.select(sp_diff_run1_n45_6Sept, C_spissum)))
sp_C_sta <- c(unlist(list.select(sp_diff_run1_n45_6Sept, C_staminodiosum)))
sp_Dicro <- c(unlist(list.select(sp_diff_run1_n45_6Sept, Dicrocaulon_sp)))
sp_Ooph <- c(unlist(list.select(sp_diff_run1_n45_6Sept, Oophytum_sp)))



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

boxplot(log(A_del), log(A_fra),  log(C_spi), log(Dicro), log(Ooph), log(D_div), log(R_com), log(R_bur), log(A_fis), log(C_sta)  )
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
