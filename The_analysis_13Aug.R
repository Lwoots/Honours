#Hi
#Final analysis script
#Finding nb variables with BRTs
#Modeling with JSDMs
#Honours 2018
#Lara Wootton

#---Set up--------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dismo, gbm, foreach, doParallel, TeachingDemos, boral, corrplot, pROC, ggplot2, dplyr)

source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

#Select study species (those that are present in ten or more plots) 
nb_sp <- species_df %>% select(plot, 
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


names(nb_sp) <- c("plot", 
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
glimpse(nb_sp)

sp_occ <- nb_sp[2:11] %>% 
  apply(2, function(x) ifelse(x == 0, 0, 1)) %>% 
  as.data.frame() %>% 
  mutate(plot = nb_sp$plot)#Create species occurrence df

#Remove unimportant soil variables and lon lat, as they aren't technically environmental

nb_soil <- my_soil_df %>% select(-c(munsell, hue, ph_h20, dN, dC, lon, lat))

#Create abundance df

abn_df <- left_join(nb_sp, nb_soil, by = "plot") %>% #Bind soil vars to species
  filter(!is.na(type)) %>% #Remove rows that we don't have data for
  select(-type, -site) #Don't need these 

abn_df <- na.omit(abn_df)

glimpse(abn_df)

#Create occurrence df

occ_df <- left_join(sp_occ, nb_soil, by = "plot") %>% #Bind soil vars to species
  filter(!is.na(type)) %>% #Remove rows that we don't have data for
  select(-type, -site) #Don't need these 

occ_site_df <- left_join(sp_occ, nb_soil, by = "plot") %>% #Bind soil vars to species
  filter(!is.na(type)) %>% #Remove rows that we don't have data for
  select(-type) #keep sites for later

occ_df <- na.omit(occ_df) #remove nas
glimpse(occ_df)

occ_site_df <- na.omit(occ_site_df)
sites <- occ_site_df$site #list of sites to filter with later

rm(sp_occ, nb_sp, nb_soil, occ_site_df)


#which variables are correlated with each other?

corr_mx <- cor(occ_df[c(12:length(occ_df))])
(big_corr <- apply(abs(corr_mx) >= 0.7 & abs(corr_mx) <1 , 1, any)) #Which correlations are greater than 0.7?
corr_vars <- occ_df %>% select(acidity,
                               Mg,
                               Na,
                               K,
                               Clay,
                               Silt,
                               Sand,
                               conductivity_ms,
                               ph_kcl, 
                               C_perc, 
                               corr_dC,
                               C_N_ratio) %>% 
  cor()

corrplot::corrplot(corr_vars, type = "lower", method = "number")

corr_vars2 <- occ_df %>% select(
                                percent_over1,
                                percent_over2,
                                Ca,
                                Na,
                                P,
                                Olsen,
                                Clay,
                                ph_kcl, 
                                N_perc,
                                corr_dN,
                                C_perc, 
                                corr_dC,
                                C_N_ratio,
                                elevation,
                                slope,
                                aspect,
                                drainage, 
                                Q_cover) %>% 
  cor()
corrplot::corrplot(corr_vars2, type = "lower", method = "number")

#Use ph as a proxy for acidity, Mg, K, corr dC (ph_et_al)
#Use Clay as a proxy for silt and sand (texture)
#Use Na as proxy for conductivity (salt)
#Use C_N_ratio as proxy for perc_C (carbon)

dat_occ <- cbind(occ_df[,1:10], occ_df %>% select(
                                                  percent_over1,
                                                  percent_over2,
                                                  Ca,
                                                  Na,
                                                  P,
                                                  Olsen,
                                                  Clay,
                                                  ph_kcl, 
                                                  N_perc,
                                                  corr_dN,
                                                  C_N_ratio,
                                                  elevation,
                                                  slope,
                                                  aspect,
                                                  drainage, 
                                                  Q_cover))
glimpse(dat_occ)
colnames(dat_occ)[colnames(dat_occ) == "ph_kcl"] <- "ph_et_al"
colnames(dat_occ)[colnames(dat_occ) == "Clay"] <- "texture"
colnames(dat_occ)[colnames(dat_occ) == "Na"] <- "salt"
colnames(dat_occ)[colnames(dat_occ) == "C_N_ratio"] <- "carbon"

#Abundance df

dat_abn <- cbind(abn_df[,2:11], abn_df %>% select(
  percent_over1,
  percent_over2,
  Ca,
  Na,
  P,
  Olsen,
  Clay,
  ph_kcl, 
  N_perc,
  corr_dN,
  C_N_ratio,
  elevation,
  slope,
  aspect,
  drainage, 
  Q_cover))

colnames(dat_abn)[colnames(dat_abn) == "ph_kcl"] <- "ph_et_al"
colnames(dat_abn)[colnames(dat_abn) == "Clay"] <- "texture"
colnames(dat_abn)[colnames(dat_abn) == "Na"] <- "salt"
colnames(dat_abn)[colnames(dat_abn) == "C_N_ratio"] <- "carbon"

glimpse(dat_abn)

rm(corr_vars, corr_vars2, big_corr, corr_mx)


#Optimise BRTs ####

#Occurrence

brt_var <- c(11:26)
response <-5
tree.com <- 10
learn <- 0.001

#Find optimal settings for BRTs - code courtesy Mike

cores <- detectCores()

#Create training function for gbm.step

step.train.fx = function(tree.com, learn){
  #set seed for reproducibility
  char2seed("StackOverflow", set = TRUE)
  k1 <- gbm.step(data = dat_abn, 
                 gbm.x = brt_var, 
                 gbm.y = response,
                 family = "poisson", 
                 tree.complexity = tree.com,
                 learning.rate = learn,
                 bag.fraction = 0.7,
                 prev.stratify=TRUE,
                 n.folds=10,
                 step.size=50,
                 silent=FALSE,
                 plot.main = TRUE,
                 n.cores= cores - 1)
  
  k.out=list(interaction.depth=k1$interaction.depth,
             shrinkage=k1$shrinkage,
             n.trees=k1$n.trees,
             AUC=k1$self.statistics$discrimination,
             cv.AUC=k1$cv.statistics$discrimination.mean,
             deviance=k1$self.statistics$mean.resid,
             cv.deviance=k1$cv.statistics$deviance.mean)  
  return(k.out)
}

#Setup to determine the optimal learning rate and tree complexity 
#define complexity and learning rate
tree.complexity<-c(2:5)
learning.rate<-c(0.01,0.005,0.001,0.0005,0.0001)

#Optimise for the current data
#setup parallel backend to use n processors
cl<-makeCluster(cores - 1)
registerDoParallel(cl)

#Run the actual function
foreach(l = tree.complexity) %do% {
  foreach(j = learning.rate) %do% {
    nam=paste0("gbm_tc",l,"lr",j)
    assign(nam,step.train.fx(tree.com=l,learn=j))
  }
}

#Stop parallel
stopCluster(cl)
registerDoSEQ()

#disable scientific notation
options(scipen=999)

#Find all item in workspace that contain "gbm_tc"
train.all<-ls(pattern="gbm_tc")

#cbind each list that contains "gbm_tc"
train.results<-list(do.call(cbind,mget(train.all)))

#Place in a data frame
train.results<- do.call(rbind, lapply(train.results, rbind))
train.results <- data.frame(matrix(unlist(train.results),ncol=7 , byrow=T))

#Change column names
colnames(train.results)<-c("TC","LR","n.trees", "AUC", "cv.AUC", "dev", "cv.dev")

#Round 4:7
train.results[,4:7]<-round(train.results[,4:7],digits=3)

#Sort by cv.dev, cv.AUC, AUC
train.results_out <-train.results[order(train.results$cv.dev,-train.results$cv.AUC, -train.results$AUC),]
names(train.results_out)

write.csv(train.results_out, "/Users/larawootton/Documents/Honours/Data/BRT_outputs/BRT_opt_Oop_occ.csv")



#Abundance





#Run BRTs ####

#Abundance
#TC = 4, LR = 0.0001

brt_var <- c(11:length(abn_occ))
tree.com <- 4
learn <- 0.0001

results <- list()
psuedoR2 <- matrix(nrow = 10, ncol = 2, dimnames = list(NULL, c("species", "R2")))


for (i in 1:10) {
  
  brt_results <- gbm.step(dat_abn,
                          gbm.x = brt_var,
                          gbm.y = i,
                          plot.main = TRUE,
                          family = "poisson",
                          step.size = 50,
                          tree.complexity = tree.com,
                          learning.rate = learn,
                          max.trees=10000,
                          n.folds = 10,
                          bag.fraction = 0.5
  )
  
  species <- names(dat_abn[i])
  
  results[[species]] <- brt_results
  
  psuedoR2[i,] <- c(species = species, R2 = 1-(brt_results$cv.statistics$deviance.mean/brt_results$self.statistics$mean.null))
  
}


results
as.data.frame(psuedoR2)

summary(results$C_spissum)
results$R_burtoniae

fin_abn_df <- dat_abn %>% 
  select(1:10, 
         ph_et_al,
         salt,
         carbon,
         Ca,
         Q_cover,
         elevation,
         percent_over1,
         percent_over2,
         P,
         N_perc,
         aspect,
         texture,
         corr_dN,
         drainage,
         slope)

#Ocurrence

#Optimising process showed that TC = 2, LR = 0.0005

brt_var <- c(11:length(dat_occ))
tree.com <- 3
learn <- 0.01


results <- list()
psuedoR2 <- matrix(nrow = 10, ncol = 2, dimnames = list(NULL, c("species", "R2")))


for (i in 1:10) {
  
  brt_results <- gbm.step(dat_occ,
                          gbm.x = brt_var,
                          gbm.y = i,
                          plot.main = TRUE,
                          family = "bernoulli",
                          step.size = 50,
                          tree.complexity = tree.com,
                          learning.rate = learn,
                          max.trees=10000,
                          n.folds = 10,
                          bag.fraction = 0.5
  )
  
  species <- names(occ_df[i])
  
  results[[species]] <- brt_results
  
  psuedoR2[i,] <- c(species = species, R2 = 1-(brt_results$cv.statistics$deviance.mean/brt_results$self.statistics$mean.null))

}


results
as.data.frame(psuedoR2)

summary(results$A_fissum)
results$R_burtoniae

#Choose variables that have 5% relative influence for at least one sp.

 fin_occ_df <- dat_occ %>% 
  select(1:10, 
         ph_et_al,
         salt,
         carbon,
         Ca,
         Q_cover,
         elevation,
         percent_over1,
         percent_over2,
         P,
         N_perc,
         aspect,
         texture,
         corr_dN,
         drainage)


 
 #JSDM ####

 #Abundance
 
 sp <- as.matrix(fin_abn_df[,1:10])
 covar <- as.matrix(fin_abn_df[,11:length(fin_abn_df)]) 
 
 abn_model_19Aug <- boral(sp,
                          X = covar,
                          family = "negative.binomial",
                          num.lv = 3,
                          save.model = T,
                          calc.ics = T)
 save(abn_model_19Aug, file = "abn_model_19Aug.rda")
 
 plot.boral(abn_model_19Aug)
 abline(0,1)
 
 #Scaled abundance
 
 sp <- as.matrix(fin_abn_df[,1:10])
 covar <- as.matrix(scale(fin_abn_df[,11:length(fin_abn_df)]))
 
 abn_model_scaled_19Aug <- boral(sp,
                          X = covar,
                          family = "negative.binomial",
                          num.lv = 3,
                          save.model = T,
                          calc.ics = T)
save(abn_model_scaled_19Aug, file = "abn_model_scaled_19Aug.rda")

plot.boral(abn_model_scaled_19Aug)
abline(0,1) 
 
summary(abn_model_scaled_19Aug)

mod_fit <- fitted.boral(abn_model_scaled_19Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)

plot(pred$R_burtoniae ~ fin_abn_df$R_burtoniae,
     ylab = "Predicted",
     xlab = "",
     main = "R_burtoniae")
abline(0,1)
plot(pred$R_comptonii ~ fin_abn_df$R_comptonii,
     xlab = "",
     ylab = "",
     main = "R_comptonii")
abline(0,1)
plot(pred$D_diversifolium ~ fin_abn_df$D_diversifolium,
     ylab = "Predicted",
     xlab = "",
     main = "D_diversifolium")
abline(0,1)
plot(pred$A_delaetii ~ fin_abn_df$A_delaetii,
     xlab = "",
     ylab = "",
     main = "A_delaetii")
abline(0,1)
plot(pred$A_fissum ~ fin_abn_df$A_fissum,
     ylab = "Predicted",
     xlab = "",
     main = "A_fissum")
abline(0,1)
plot(pred$A_framesii ~ fin_abn_df$A_framesii,
     xlab = "",
     ylab = "",
     main = "A_framesii")
abline(0,1)
plot(pred$C_spissum ~ fin_abn_df$C_spissum,
     ylab = "Predicted",
     xlab = "",
     main = "C_spissum")
abline(0,1)
plot(pred$C_staminodiosum ~ fin_abn_df$C_staminodiosum,
     xlab = "",
     ylab = "",
     main = "C_staminodiosum")
abline(0,1)
plot(pred$Dicrocaulon_sp ~ fin_abn_df$Dicrocaulon_sp,
     ylab = "Predicted",
     xlab = "Observed",
     main = "Dicrocaulon_sp")
abline(0,1)
plot(pred$Oophytum_sp ~ fin_abn_df$Oophytum_sp,
     ylab = "Observed",
     xlab = "",
     main = "Oophytum_sp")
abline(0,1)


#Full model occurrence
sp <- as.matrix(fin_occ_df[,1:10])
covar <- as.matrix(fin_occ_df[,11:length(fin_occ_df)]) 
 
#occ_model_14Aug <- boral(sp,
#                         X = covar,
#                         family = "binomial",
#                         num.lv = 3,
#                         save.model = T,
#                         calc.ics = T)

save(occ_model_14Aug, file = "occ_model_14Aug.rda")
load("occ_model_14Aug.rda")

plot.boral(occ_model_14Aug)
summary(occ_model_14Aug)
lvsplot(occ_model_14Aug)
coefsplot("ph_et_al", occ_model_14Aug)

var <- calc.varpart(occ_model_14Aug)
part <- data.frame(enviro = var$varpart.X, bio = var$varpart.lv)

barplot(as.matrix(t(part)), las = 2, col = c("brown4", "light green"))


mod_fit <- fitted.boral(occ_model_14Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)
pred_all <- ROCR::prediction(pred, fin_occ_df[,1:10])
rocs <- ROCR::performance(pred_all, "tpr", "fpr")

colours <- c("burlywood", "firebrick2", "chocolate1","darkolivegreen1", "darkseagreen1", "green4", "orange4", "yellow2", "darkblue","turquoise4"  )
plot(rocs, col = as.list(colours))
abline(0,1)

aucs <- ROCR::performance(pred_all, "auc")

#scaled vars to improve model fit

sp <- as.matrix(fin_occ_df[,1:10])
covar <- as.matrix(scale(fin_occ_df[,11:24])) 

occ_model_scaled_nospat_17Aug <- boral(sp,
                         X = covar,
                         family = "binomial",
                         num.lv = 3,
                         save.model = T, 
                         calc.ics = T)

save(occ_model_scaled_nospat_17Aug, file = "occ_model_scaled_nospat_17Aug.rda")


plot.boral(occ_model_scaled_nospat_17Aug)
summary(occ_model_scaled_nospat_17Aug)
lvsplot(occ_model_scaled_nospat_17Aug)


par(mfrow = c(2,3))
coefsplot("percent_over1", occ_model_scaled_nospat_17Aug)
coefsplot("ph_et_al", occ_model_scaled_14Aug)
coefsplot("percent_over1", occ_model_scaled_14Aug)
coefsplot("carbon", occ_model_scaled_14Aug)
coefsplot("texture", occ_model_scaled_14Aug)
coefsplot("Ca", occ_model_scaled_14Aug)


var <- calc.varpart(occ_model_scaled_nospat_17Aug)
part <- data.frame(enviro = var$varpart.X, bio = var$varpart.lv)
barplot(as.matrix(t(part)), las = 2, col = c("brown4", "light green"))

bardat <- data.frame(Species =rep(rownames(part),2), variance = c(as.vector(part[,1]), as.vector(part[,2])), type = c(rep("enviro",10), rep("lv",10)))
glimpse(bardat)

ggplot(bardat, aes(x=Species, y=variance, fill = type)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   colour = "black", 
                                   size = 10,
                                   vjust = 0.98, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  labs(y = "Proportion of variance")
  

?element_text

envcors <- get.enviro.cor(occ_model_scaled_nospat_17Aug)
rescors <- get.residual.cor(occ_model_scaled_nospat_17Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor, order = "hclust")
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)


qgraph::qgraph(envcors$sig.cor, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10)
qgraph::qgraph(rescors$sig.correlaton, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10)

#Caterpillar plots


mod <- occ_model_scaled_nospat_17Aug
par(mfrow = c(1,2))

for (i in 1:10) {
  col.seq <-
    rep("black", length(mod$hpdintervals$X.coefs[i, 1:14, "lower"]))
  col.seq[mod$hpdintervals$X.coefs[i, 1:14, "lower"] < 0 &
            mod$hpdintervals$X.coefs[i, 1:14, "upper"] > 0] <- "grey"
  
  plot(
    x = c(mod$X.coefs.mean[i, 1:14]),
    y = 1:14,
    yaxt = "n",
    ylab = "",
    xlab = rownames(mod$X.coefs.mean)[i],
    col = col.seq,
    xlim = c(
      min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
      max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
    ),
    pch = "x"
  )
  
  axis(
    2,
    labels = colnames(mod$X.coefs.mean),
    at = 1:14,
    las = 2,
    cex.axis = 0.75
  )
  
  segments(
    x0 = mod$hpdintervals$X.coefs[i, 1:14, "lower"],
    y0 = 1:14,
    x1 = mod$hpdintervals$X.coefs[i, 1:14, "upper"],
    y1 = 1:14,
    col = col.seq
  )
  abline(v = 0, lty = 3)
  
}


#Predicting occurrence with model ####

mod_fit <- fitted.boral(occ_model_scaled_nospat_17Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)
summary(pred)
plot(pred$C_spissum ~ jitter(fin_occ_df$C_spissum))
roc1 <- roc(fin_occ_df$C_spissum, pred$C_spissum)
plot(roc1)

#Predicting occurrence onto new data ####
#Random 50% of plots
rplots_19aug <- sample(150, 100)
#rplots <- c(76, 37, 32,28,116,23, 48,16, 146, 45,18, 132,88, 5, 43, 126, 59, 118,  38,  31,  93, 123, 136, 100,3,57,  86,80,  13, 115, 112, 1, 33, 150,  68, 72,96, 41, 95,  15,  63, 7, 29,  17,20,40, 109, 142,51, 67,19,79, 71,56, 149, 139, 65, 39, 6, 114, 104, 141, 145, 92, 49,55, 46, 130,66, 144, 60,61, 137,  47,  87)
save(rplots_19aug, file = "rplots_19aug.rda")

train_occ <- fin_occ_df[rplots_19aug,]

sp <- as.matrix(train_occ[,1:10])
covar <- as.matrix(scale(train_occ[,11:24])) 

train_occ_scaled_19Aug <- boral(sp,
                                X = covar,
                                family = "binomial",
                                num.lv = 3,
                                save.model = T)

save(train_occ_scaled_19Aug, file = "train_occ_scaled_19Aug.rda")
load("train_occ_scaled_19Aug.rda")
plot.boral(train_occ_scaled_19Aug)

mod_fit <- fitted.boral(train_occ_scaled_19Aug, est = "mean")
pred <- as.data.frame(mod_fit$out)
summary(pred)
plot(pred$C_spissum ~ jitter(train_occ$C_spissum))
roc1 <- roc(train_occ$C_spissum, pred$C_spissum)
plot(roc1)

test_occ <- fin_occ_df[-rplots_19aug,]
test_covar <- as.matrix(scale(test_occ[,11:24])) 

newpred <- predict.boral(train_occ_scaled_19Aug, 
                         newX = test_covar, 
                         predict.type = "marginal",
                         est = "mean")

newpred$linpred
newpred$lower

pred <- newpred$linpred
pred_all <- ROCR::prediction(pred, test_occ[,1:10])
rocs <- ROCR::performance(pred_all, "tpr", "fpr")
plot(rocs, col = as.list(colours))
abline(0,1)

aucs <- ROCR::performance(pred_all, "auc")


#Use two plots to predict onto one

#2 and 3 to 1
site_occ <- cbind(site = sites, fin_occ_df)
plots_23 <- site_occ %>% 
  filter(!site == "site1")
plots_1 <- site_occ %>% 
  filter(site == "site1")


sp <- as.matrix(plots_23[,2:11])
covar <- as.matrix(scale(plots_23[,12:25])) 

plots23_19Aug <- boral(sp,
                                X = covar,
                                family = "binomial",
                                num.lv = 3,
                                save.model = T)
save(plots23_19Aug, file = "plots23_19Aug.rda")

plot.boral(plots23_15Aug)

newpred <- predict.boral(plots23_19Aug, 
                         newX = plots_1[,12:25], 
                         predict.type = "marginal",
                         est = "mean")
summary(newpred$linpred)

pred <- ROCR::prediction(newpred$linpred[,4], plots_1[,5])
rocs <- ROCR::performance(pred, "tpr", "fpr")
(auc <- ROCR::performance(pred, "auc"))
plot(rocs)

preds_list <- list(newpred$linpred[,1], newpred$linpred[,2], newpred$linpred[,3],newpred$linpred[,4],newpred$linpred[,5], newpred$linpred[,6],newpred$linpred[,7],newpred$linpred[,8],newpred$linpred[,9],newpred$linpred[,10])

actuals_list <- list(plots_1[,2],plots_1[,3],plots_1[,4],plots_1[,5],plots_1[,6],plots_1[,7],plots_1[,8],plots_1[,9],plots_1[,10],
                     plots_1[,11])

pred <- ROCR::prediction(preds_list, actuals_list)
rocs <- ROCR::performance(pred, "tpr", "fpr")
aucs <- ROCR::performance(pred, "auc")
plot(rocs, col = as.list(1:10))

#1 and 2 to 3

plots_12 <- site_occ %>% 
  filter(!site == "site3")
plots_3 <- site_occ %>% 
  filter(site == "site3")

sp <- as.matrix(plots_12[,2:11])
covar <- as.matrix(scale(plots_12[,12:25])) 

plots12_19Aug <- boral(sp,
                       X = covar,
                       family = "binomial",
                       num.lv = 3,
                       save.model = T)
save(plots12_19Aug, file = "plots12_19Aug.rda")

plot.boral(plots12_15Aug)

newpred2 <- predict.boral(plots12_17Aug, 
                         newX = plots_3[,12:25], 
                         predict.type = "marginal",
                         est = "mean")

pred <- newpred2$linpred

pred_all <- ROCR::prediction(pred[,c(1:5, 7:10)], plots_3[,c(2:6,8:11)])
rocs <- ROCR::performance(pred_all, "tpr", "fpr")
plot(rocs, col = as.list(colours))
abline(0,1)

aucs <- ROCR::performance(pred_all, "auc")

#sites 1 and 3 on 2
plots_13 <- site_occ %>% 
  filter(!site == "site2")
plots_2 <- site_occ %>% 
  filter(site == "site2")

sp <- as.matrix(plots_13[,2:11])
covar <- as.matrix(scale(plots_13[,12:25])) 

plots13_19Aug <- boral(sp,
                       X = covar,
                       family = "binomial",
                       num.lv = 3,
                       save.model = T)

save(plots13_19Aug, file = "plots13_1Aug.rda")

plot.boral(plots13_15Aug)

newpred3 <- predict.boral(plots13_19Aug, 
                          newX = plots_2[,12:25], 
                          predict.type = "marginal",
                          est = "mean")

pred <- newpred3$linpred

pred_all <- ROCR::prediction(pred[,c(1, 3:7, 9, 10)], plots_2[,c(2,4:8, 10,11)])
rocs <- ROCR::performance(pred_all, "tpr", "fpr")
plot(rocs, col = as.list(colours))
abline(0,1)

aucs <- ROCR::performance(pred_all, "auc")


#Predicting abundance
plots_2[2:11]



