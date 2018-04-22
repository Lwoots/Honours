#First attempt at kriging
#Honours 2018
#22 April
#Resources: https://rpubs.com/nabilabd/118172, http://rstudio-pubs-static.s3.amazonaws.com/80464_9156596afb2e4dcda53e3650a68df82a.html

#------------------------------------------------------------
rm(list = ls())
setwd("/Users/larawootton/Documents/Honours/Data")

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, RColorBrewer, sp, gstat, automap)

ph <- read.csv("pH_for_analysis.csv", sep = ";")
raw <- read.csv("Species_abundance_for_analysis.csv", sep = ";")

#Data wrangling ####

#Extract the names of all the plots
plots <- data.frame(raw[,2])
names(plots) <- "plot"
rm(raw)

#Add pseudo coordinates and site names to the full dataset
full_grid <-  plots %>% mutate(y = as.vector(rep(sort(c(rep(1,11), rep(2,11), rep(3,11), rep(4,11), rep(5,11), rep(6,11), rep(7,11), rep(8,11), rep(9,11), rep(10,11), rep(11,11)), decreasing = T), 3)),
                             x = as.vector(rep(rep(seq(1, 11, 1), 11), 3)), site = c(rep("site1", 121), rep("site2", 121), rep("site3",121)))

glimpse(full_grid)
plot(full_grid$x, full_grid$y)

#now get the smaller grid by joining ph to the full dataset

measured <- left_join(ph, full_grid, by = "plot")
glimpse(measured)

#Extact the full grid of plots from site1
site1_full_grid <- full_grid %>% filter(site == "site1")

#Extract small grid from site 1

site1_measured <- measured %>% filter(site == "site1", type == "grid")



#Now for the actual kriging####

#convert to spatial dataframe

coordinates(site1_measured) <- c("x", "y")

#Create a semivariogram

ph_vgm <- variogram(ph_kcl ~ 1, site1_measured)
plot(ph_vgm) #From here you could fit your own model using fit.variogram

#Fitting a model automatically

auto_ph <- autofitVariogram(ph_kcl ~ 1, site1_measured)
summary(auto_ph)

plot(auto_ph)

#To auto Krige
coordinates(site1_full_grid) <- c("x", "y")

ph_autoKrige <- autoKrige(ph_kcl ~ 1, site1_measured, site1_full_grid)

plot(ph_autoKrige)
glimpse(ph_autoKrige)
ph_autoKrige@var1.pred

pred_krige <- as.data.frame(ph_autoKrige[1])

p <- ph_autoKrige[1] %>% as.data.frame %>%
  ggplot(aes(x=krige_output.x, y=krige_output.y)) + geom_tile(aes(fill=krige_output.var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  #scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() +
  theme(legend.position = "none")

test <- measured %>% filter(site == "site1", type == "random")

p2 <- ggplot(test,
       aes(x,y)
       ) + 
  geom_tile(aes(fill = ph_kcl)) +
  scale_fill_gradient(low = "yellow", high="red", limits = c(min(pred_krige$krige_output.var1.pred), max(pred_krige$krige_output.var1.pred)))

t.test(test$ph_kcl, pred_krige$krige_output.var1.pred)
plot_grid(p,p2)


#Try with logging

log_auto_ph <- autofitVariogram(log(ph_kcl) ~ 1, site1_measured)
log_ph_autoKrige <- autoKrige(log(ph_kcl) ~ 1, site1_measured, site1_full_grid)
log_pred_krige <- as.data.frame(log_ph_autoKrige[1])

p <- log_ph_autoKrige[1] %>% as.data.frame %>%
  ggplot(aes(x=krige_output.x, y=krige_output.y)) + geom_tile(aes(fill=krige_output.var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  #scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw() +
  theme(legend.position = "none")

p2 <- ggplot(test,
             aes(x,y)
) + 
  geom_tile(aes(fill = log(ph_kcl))) +
  scale_fill_gradient(low = "yellow", high="red", limits = c(min(log_pred_krige$krige_output.var1.pred), max(log_pred_krige$krige_output.var1.pred)))

t.test(test$ph_kcl, pred_krige$krige_output.var1.pred)
plot_grid(p,p2)
