#Final analysis script
#Finding nb variables with BRTs
#Modeling with JSDMs
#Honours 2018
#Lara Wootton

#---Set up--------------------------------
rm(list = ls())

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, dismo, gbm, foreach, doParallel, TeachingDemos)

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

#Remove unimportant soil variables

nb_soil <- my_soil_df %>% select(-c(munsell, hue, ph_h20, dN, dC))

#Create abundance df

abn_df <- left_join(nb_sp, nb_soil, by = "plot") %>% #Bind soil vars to species
  filter(!is.na(type)) %>% #Remove rows that we don't have data for
  select(-type, -site) #Don't need these 

glimpse(abn_df)

#Create occurrence df

occ_df <- left_join(sp_occ, nb_soil, by = "plot") %>% #Bind soil vars to species
  filter(!is.na(type)) %>% #Remove rows that we don't have data for
  select(-type, -site) #Don't need these 

rm(sp_occ, nb_sp, nb_soil)

#Optimise BRTs ####

#Occurrence

brt_var <- c(12:45)
response <- 4
tree.com <- 10
learn <- 0.001

#Find optimal settings for BRTs - code courtesy Mike

cores <- detectCores()

#Create training function for gbm.step

step.train.fx = function(tree.com, learn){
  #set seed for reproducibility
  char2seed("StackOverflow", set = TRUE)
  k1 <- gbm.step(data = occ_df, 
                 gbm.x = brt_var, 
                 gbm.y = response,
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

#Ocurrence


#Abundance



#JSDM ####


#Predicting occurrence with model ####

#Predicting occurrence onto new data ####

#Predicting abundance



