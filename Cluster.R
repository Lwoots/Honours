#Using dice for cluster analysis
#Are there different communities?
#Honours 2018
#10 June

#---------------------------------------------------------------------------

rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, vegan, ape)

source("/Users/larawootton/Documents/Honours/Project_analysis/to_be_sourced.R")

nb_sp <- species_df %>% select(plot, ruschia_burtoniae, c_spissum, a_delaetii, oophytum, a_fissum, dicrocaulon, drosanthemum_diversifolium, a_framesii, mesemb_1, c_subfenestratum, co_calculus, galenia_fruticosa, crassula_muscosa, tylecodon_pygmaeus, cephalophyllum_staminodiosum)

row.names(nb_sp) <- nb_sp$plot

nb_sp <- nb_sp %>% select(-plot)

sp_occ <- nb_sp %>% apply(2, function(x) ifelse(is.na(x), 0, 1))

dist <- designdist(sp_occ,method="(A+B-2*J)/(A+B)",terms="binary",abcd=FALSE)
hclust(dist,method="average")

plot(dist)
