rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")
source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, dismo, gbm, foreach, doParallel, TeachingDemos)

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


#Find optimal settings for BRTs - code courtesy Mike

cores <- detectCores()

#Create training function for gbm.step

step.train.fx = function(tree.com, learn){
  #set seed for reproducibility
  char2seed("StackOverflow", set = TRUE)
  k1 <- gbm.step(data = my_prez, 
                 gbm.x = brt_var, 
                 gbm.y = 7,
                 family = "bernoulli", 
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
tree.complexity<-c(1:5)
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

write.csv(train.results_out, "/Users/larawootton/Documents/Honours/Data/BRT_outputs/A_delaetii_ben_train_results_out")