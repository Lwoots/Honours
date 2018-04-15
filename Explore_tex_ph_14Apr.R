#Honours 2018
#Exploratory analysis 14 April
#Descriptive stats for texture and ph, conductivity

#----------------------------------------------------#
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, visdat, ggplot2, RColorBrewer, cowplot)

plots <- read.csv("plots.csv")
sieved <- read.csv("Soil_texture_Final.csv", sep = ";") 
malvern <- read.csv("Malvern_dodge_labels.csv", sep = ";")
ph <- read.csv("ph_cond_13Apr.csv", sep = ";")


which(duplicated(malvern$Sample) == TRUE) #Find the duplicated samples in the malvern data

setdiff(plots$plot, sieved$plot)
setdiff(plots$plot, malvern$Sample)
setdiff(malvern$Sample, sieved$plot)

names(malvern)[1] <- "plot"
names(ph)[1] <- "plot"

texture <- left_join(sieved, malvern, by = "plot")

texture <- texture %>% select(1:12, 127:133)
glimpse(texture)

texture <- left_join(texture, ph, by = "plot")
texture <- mutate(texture, 
                  Sand = Fine_sand + Medium_sand + Coarse_sand + Very_coarse_sand,
                  prop_less2 = (between2_1/soil_mass)*100, 
                  prop_less5 = (between2_5mm/soil_mass)*100)

grid_texture <- texture %>% filter(type == "grid") %>% arrange(plot)

site3 <- grid_texture[7:42,] %>% 
                     mutate(y = sort(as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6))), decreasing = T),
                            x = as.vector(rep(seq(1, 11, 2), 6)))

site1 <- grid_texture[c(1:6, 43:72),] %>% 
                      mutate(y = sort(as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6))), decreasing = T),
                             x = as.vector(rep(seq(1, 11, 2), 6)))

site2 <- grid_texture[73:111,]

rm(ph, malvern, plots, sieved)

#Descriptive stats ####
names(texture)
clay <- texture %>% summarise(avg = mean(Clay), 
                              max = max(Clay), 
                              min = min(Clay), 
                              sd = sd(Clay))
silt <- texture %>% summarise(avg = mean(Silt), 
                              max = max(Silt), 
                              min = min(Silt), 
                              sd = sd(Silt))   
very_fine_sand <- texture %>% summarise(avg = mean(Very_fine_sand), 
                                        max = max(Very_fine_sand), 
                                        min = min(Very_fine_sand), 
                                        sd = sd(Very_fine_sand))
medium_sand <- texture %>% summarise(avg = mean(Medium_sand), 
                                     max = max(Medium_sand), 
                                     min = min(Medium_sand), 
                                     sd = sd(Medium_sand))
coarse_sand <- texture %>% summarise(avg = mean(Coarse_sand), 
                                     max = max(Coarse_sand), 
                                     min = min(Coarse_sand), 
                                     sd = sd(Coarse_sand))
very_coarse_sand <- texture %>% summarise(avg = mean(Very_coarse_sand), 
                                          max = max(Very_coarse_sand), 
                                          min = min(Very_coarse_sand), 
                                          sd = sd(Very_coarse_sand))
total_sand <- texture %>% summarise(avg = mean(Sand), 
                                 max = max(Sand), 
                                 min = min(Sand), 
                                 sd = sd(Sand))
less2 <- texture %>% summarise(avg = mean(prop_less2), 
                               max = max(prop_less2), 
                               min = min(prop_less2), 
                               sd = sd(prop_less2))
less5 <- texture %>% summarise(avg = mean(prop_less5),
                               
                               max = max(prop_less5), min = min(prop_less5), 
                               sd = sd(prop_less5))

tex_des <- rbind(clay, silt, very_fine_sand, medium_sand, coarse_sand, very_coarse_sand, total_sand, less2, less5)
category <- c("clay", "silt", "very_fine_sand", "medium_sand", "coarse_sand", "very_coarse_sand","total_sand", "less2", "less5")
data.frame(category, tex_des)

rm(clay, silt, very_fine_sand, medium_sand, coarse_sand, very_coarse_sand, total_sand, less2, less5)



#Plots####
theme_set(theme_bw())
ggplot(site3,
       aes(x,y, colour = Clay)) +
  scale_color_gradient(low = "blue", high = "red", limits = c(0,8))+
  geom_point()
  


my_theme <- theme(
  #axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid = element_blank()
)
my_colour <- scale_color_distiller(palette = "YlOrRd", limits = c(0,8.03), direction = 0)
#Clay
names(site1)
leg <- ggplot(site1,
              aes(x, y,
                  colour = Clay)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1")
p1 <- ggplot(site1,
             aes(x, y,
                 colour = Clay)) +
              geom_point(cex = 9) +
              my_theme +
              my_colour +
              labs(x = "Site 1") +
  theme(legend.position = "none")
p2 <- ggplot(site3,
             aes(x, y,
                 colour = Clay)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 3") +
  theme(legend.position = "none")


legend <- get_legend(leg)
plot_grid(p1,p2, legend, nrow = 1, rel_widths = c(1,1,0.1))            



#Silt
my_colour <- scale_color_distiller(palette = "YlOrRd", limits = c(min(grid_texture$Silt), max(grid_texture$Silt)), direction = 0)

leg <- ggplot(site1,
              aes(x, y,
                  colour = Silt)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1")
p1 <- ggplot(site1,
             aes(x, y,
                 colour = Silt)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1") +
  theme(legend.position = "none")
p2 <- ggplot(site3,
             aes(x, y,
                 colour = Silt)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 3") +
  theme(legend.position = "none")


legend <- get_legend(leg)
plot_grid(p1,p2, legend, nrow = 1, rel_widths = c(1,1,0.1))   

#Sand

my_colour <- scale_color_distiller(palette = "YlOrRd", limits = c(min(grid_texture$Sand), max(grid_texture$Sand)), direction = 0)

leg <- ggplot(site1,
              aes(x, y,
                  colour = Sand)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1")
p1 <- ggplot(site1,
             aes(x, y,
                 colour = Sand)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1") +
  theme(legend.position = "none")
p2 <- ggplot(site3,
             aes(x, y,
                 colour = Sand)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 3") +
  theme(legend.position = "none")


legend <- get_legend(leg)
plot_grid(p1,p2, legend, nrow = 1, rel_widths = c(1,1,0.1)) 

#ph

my_colour <- scale_color_distiller(palette = "YlOrRd", limits = c(min(grid_texture$ph_kcl), max(grid_texture$ph_kcl)), direction = 0)

leg <- ggplot(site1,
              aes(x, y,
                  colour = ph_kcl)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1")
p1 <- ggplot(site1,
             aes(x, y,
                 colour = ph_kcl)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 1") +
  theme(legend.position = "none")
p2 <- ggplot(site3,
             aes(x, y,
                 colour = ph_kcl)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 3") +
  theme(legend.position = "none")

my_colour <- scale_color_distiller(palette = "YlOrRd", limits = c(min(grid_texture$conductivity_ms), max(grid_texture$conductivity_ms)), direction = 0)


ggplot(site3,
       aes(x, y,
           colour = conductivity_ms)) +
  geom_point(cex = 9) +
  my_theme +
  my_colour +
  labs(x = "Site 3")


plot(site3$x, site3$y)
text(site3$x, site3$y, labels = site3$plot)
plot(site3$y, site3$x)
text(site3$y, site3$x, labels = site3$plot)

?text
