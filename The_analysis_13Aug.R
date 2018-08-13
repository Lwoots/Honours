#Final analysis script
#Finding nb variables with BRTs
#Modeling with JSDMs
#Honours 2018
#Lara Wootton

#---Set up--------------------------------
rm(list = ls())
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

#Abundance



#Run BRTs ####

#Ocurrence


#Abundance



#JSDM ####


#Predicting occurrence with model ####

#Predicting occurrence onto new data ####

#Predicting abundance




