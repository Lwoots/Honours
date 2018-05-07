#Honours 2018
#Exploratory analysis 19 April
#Spatial distribution of abundance

#------------------------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, visdat, ggplot2, RColorBrewer, cowplot)

raw <- read.csv("Species_abundance_for_analysis.csv", sep = ";")

count <- function(x, na.rm = F) {
  ifelse(x > 0, 1, 0)
}
raw <- raw %>% mutate(site = c(rep("site1", 121), rep("site2", 121), rep("site3",121)))




vis_dat(raw)
raw$site <- as.factor(raw$site)
summary(raw$plot)
str(raw$ruschia_burtoniae)

plot(raw$site, raw$monilaria_moniliformis)
plot(raw$site, raw$c_spissum)
plot(raw$co_calculus ~ raw$site)
plot(raw$a_fissum ~ raw$site)
plot(raw$c_subfenestratum ~ raw$site)
plot(raw$oophytum ~ raw$site)
plot(raw$dicrocaulon ~ raw$site)
plot(raw$drosanthemum_diversifolium ~ raw$site)
plot(raw$crassula_deceptor ~ raw$site)
plot(raw$mesemb_1 ~ raw$site)
plot(raw$tylecodon_pygmaeus ~ raw$site)
plot(raw$Crassula.sp..mini.columnus ~ raw$site)
plot(raw$monoleria_pisiformis ~ raw$site)
plot(raw$framesiiXD ~ raw$site)
plot(raw$crassula_muscosa ~ raw$site)
plot(raw$sarcocornia ~ raw$site)
plot(raw$sarcocaulon_spikey ~ raw$site)
plot(raw$caroxylon_sp1 ~ raw$site)
plot(raw$cephalophyllum_staminodiosum ~ raw$site)
plot(raw$Crassula.sp..mini.columnus ~ raw$site)


grid <- raw %>% mutate(y = as.vector(rep(sort(c(rep(1,11), rep(2,11), rep(3,11), rep(4,11), rep(5,11), rep(6,11), rep(7,11), rep(8,11), rep(9,11), rep(10,11), rep(11,11)), decreasing = T), 3)),
                       x = as.vector(rep(rep(seq(1, 11, 1), 11), 3)))

plot(grid$x, grid$y)

text(grid$x[grid$site == "site2"], grid$y[grid$site == "site2"], labels = grid$plot[grid$site == "site3"])


site1 <- grid %>% filter(site == "site1")
site2 <- grid %>% filter(site == "site2")
site3 <- grid %>% filter(site == "site3")


#Theme ####
theme_set(theme_bw())
my_theme <- theme(
  #axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid = element_blank(),
  legend.title.align = 0.5,
  legend.title = element_text(angle = 270, hjust = 1)
)

#Plots####

my_colour <- scale_color_distiller(palette = "YlGn", limits = c(min(grid$ruschia_burtoniae, na.rm = T), max(grid$ruschia_burtoniae, na.rm = T)), direction = 0, name = "ruschia_burtoniae")

leg <- ggplot(site1,
              aes(x, y,
                  colour = ruschia_burtoniae)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour +
  labs(x = "Site 1")
p1 <- ggplot(site1,
             aes(x, y,
                 colour = ruschia_burtoniae)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour +
  labs(x = "Site 1") +
  theme(legend.position = "none") 

p2 <- ggplot(site2,
             aes(x, y,
                 colour = ruschia_burtoniae)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour +
  labs(x = "Site 2") +
  theme(legend.position = "none")

p3 <- ggplot(site3,
             aes(x, y,
                 colour = ruschia_burtoniae)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour +
  labs(x = "Site 3") +
  theme(legend.position = "none")


legend <- get_legend(leg)
plot_grid(p1,p2,p3, legend, nrow = 1, rel_widths = c(0.9,0.9,0.9,0.2)) 


#Species richness ####

raw2 <- read.csv("Species_abundance_for_analysis.csv", sep = ";")
species_only <- raw2 %>% select(-c(Q_cover, transect, X, lampranthus)) 
glimpse(species_only)

species_richness<- replace(species_only[,2:46], species_only[,2:46] > 0, 1)
glimpse(species_richness)

species_richness <- apply(species_richness, 1, sum, na.rm = T )


sp_richness <- data.frame(plot = species_only[,1], sp_rich = species_richness, site = c(rep("site1", 121), rep("site2", 121), rep("site3",121)))
glimpse(sp_richness)

sp_richness <- sp_richness %>% mutate(y = as.vector(rep(sort(c(rep(1,11), rep(2,11), rep(3,11), rep(4,11), rep(5,11), rep(6,11), rep(7,11), rep(8,11), rep(9,11), rep(10,11), rep(11,11)), decreasing = T), 3)),
                                              x = as.vector(rep(rep(seq(1, 11, 1), 11), 3)))
glimpse(sp_richness)


my_colour2 <- scale_color_gradient(low = "white", high = "red", limits = c(min(sp_richness$sp_rich, na.rm = T), max(sp_richness$sp_rich, na.rm = T)), breaks = seq(0,9,1), name = "Species richness")

leg <- ggplot(sp_richness[sp_richness$site == "site1",],
              aes(x, y,
                  colour = sp_rich)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour2 +
  labs(x = "Site 1")
p1 <- ggplot(sp_richness[sp_richness$site == "site1",],
             aes(x, y,
                 colour = sp_rich)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour2 +
  labs(x = "Site 1") +
  theme(legend.position = "none") 

p2 <- ggplot(sp_richness[sp_richness$site == "site2",],
             aes(x, y,
                 colour = sp_rich)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour2 +
  labs(x = "Site 2") +
  theme(legend.position = "none")

p3 <- ggplot(sp_richness[sp_richness$site == "site3",],
             aes(x, y,
                 colour = sp_rich)) +
  geom_point(cex = 10, shape = 15) +
  my_theme +
  my_colour2 +
  labs(x = "Site 3") +
  theme(legend.position = "none")


legend <- get_legend(leg)
plot_grid(p1,p2,p3, legend, nrow = 1, rel_widths = c(0.9,0.9,0.9,0.2)) 
