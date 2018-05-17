#Script to tie together all the different datasets

#Honours 2018

setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, munsell)

plots <- read.csv("geo_ref_plots_15May.csv")
sieved <- read.csv("Soil_texture_Final.csv", sep = ";") 
ph <- read.csv("pH_for_analysis.csv", sep = ";")
raw_els <- read.csv("Elsenburg_data_for_analysis.csv", sep = ";")
munsell <- read.csv("Munsell_all_plots_for_analysis.csv", sep = ";")
texture <- read.csv("Malvern_dodge_labels.csv", sep = ";")
abundance <- read.csv("Species_abundance_for_analysis.csv", sep = ";")


#Species abundance dataframe

sp_abundance <- left_join(plots, abundance,by = "plot")
species_only <- sp_abundance %>% select(-c(Q_cover, transect, X, lampranthus, site, lon, lat)) 
species_richness<- replace(species_only[,2:46], species_only[,2:46] > 0, 1)
species_richness <- apply(species_richness, 1, sum, na.rm = T )

species_df <- data.frame(sp_abundance, species_richness) 

#getting rid of typos

species_df$sarcocaulon_spikey[258] <- 2
species_df$framesiiXD[234] <- NA
species_df$a_framesii[234] <- 7

rm(species_only, sp_abundance, abundance, species_richness)


#Soil data_frame
levels(raw_els$plot)[levels(raw_els$plot)== "R40"] <- "R40 "
levels(raw_els$plot)[levels(raw_els$plot)== "R0"] <- "R0 "
levels(raw_els$plot)[levels(raw_els$plot)== "R10"] <- "R10 "
levels(raw_els$plot)[levels(raw_els$plot)== "R20"] <- "R20 "
levels(raw_els$plot)[levels(raw_els$plot)== "R30"] <- "R30 "
levels(raw_els$plot)[levels(raw_els$plot)== "R50"] <- "R50 "

texture <- texture %>% select(-c(2:113, 121))
texture <- mutate(texture, 
                  Sand = Fine_sand + Medium_sand + Coarse_sand + Very_coarse_sand)
raw <- left_join(plots, sieved, by = "plot")

raw <- raw %>% select(-c(5:13))

munsell_150 <- sieved %>% select(plot, munsell)
munsell_150$munsell[munsell_150$munsell == "10YR 7/3"] <- "10YR 7/4"
munsell_150$munsell[munsell_150$munsell == "7.5YR 6/3"] <- "7.5YR 6/4"
munsell_remainder <- read.csv("Munsell_all_plots_for_analysis.csv", sep = ";", na.strings = "")



mun <- munsell_remainder[,c(1,4)]
names(mun) <- c("plot", "munsell")
tog <- rbind(mun, munsell_150)
full_munsell <- left_join(raw, tog, by = "plot")
full_munsell <- left_join(full_munsell, munsell_remainder, by = "plot")
full_munsell <- full_munsell %>% select(-c(munsell.x, date, munsell, even))
names(full_munsell)[8] <- "munsell"

break_down <- mnsl2hvc(full_munsell$munsell)
colour <- data.frame(full_munsell, break_down)
full_dat <- colour %>% mutate(my_hue = if_else(hue == "5YR", 5, if_else(hue == "7.5YR", 7.5, if_else(hue == "10YR", 10, 0))))

more_dat <- left_join(full_dat, raw_els, by = "plot")



names(texture)[1] <- "plot"
tex <- left_join(more_dat, texture, by = "plot")

my_soil_df <- tex %>% mutate(Sand = Fine_sand + Medium_sand + Coarse_sand + Very_coarse_sand)
my_soil_df <- left_join(my_soil_df, ph[,-c(2,6,7)], by = "plot")

my_soil_df$type[c(178, 180, 200, 202)] <- "random"
my_soil_df$type[190] <- "grid"

rm(break_down, colour, full_dat, full_munsell, more_dat, mun, munsell_150, munsell_remainder, ph, plots, raw, raw_els, sieved, tex, texture, tog, munsell)