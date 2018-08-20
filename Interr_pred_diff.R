#Interrogating prediction differences
#Cause of discrepancy between predictive ability
#Honours 2018

#Set up -----------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(ggplot2, MASS, dplyr)

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
site_sums <- rbind(site_sums, tot = total)

sum_dat <- as.data.frame(t(site_sums[,2:11]))
colnames(sum_dat) <- c("Site 1", "Site 2", "Site 3", "Total")
sum_dat