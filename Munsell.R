
if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(dplyr, ggplot2, munsell)

setwd("/Users/larawootton/Documents/Honours/Data")

### Creating the soil colour grids 21 March ###


soil_col <- read.csv("Soil_texture_Final.csv", sep = ";")
grid <- subset(soil_col, type == "grid")

rownames(grid) <- 1:length(grid$munsell)

plot_mnsl(grid$munsell) #Doesn't like uneven numbers

grid$munsell[grid$munsell == "10YR 7/3"] <- "10YR 7/4"
grid$munsell[grid$munsell == "7.5YR 6/3"] <- "7.5YR 6/4"



#Site 3####

plot_mnsl(grid$munsell[c(15:36,2,37, 3:14)])

#Site 1 ####

plot_mnsl(grid$munsell[38:73])

#Site 2 ####
dat <- data.frame(rbind(grid$munsell[c(74:90,103:106)], NA, grid$munsell[c(107, 91:102)]))

temp <- matrix(rep(NA,14), nrow = 1)
colnames(temp) <- colnames(grid)
dat <- rbind(grid[c(74:84, 1, 85:90,103:106),],temp, grid[c(107, 91:102),])

plot_mnsl(dat$munsell)


plot_mnsl(grid$munsell[1:20], na.rm = T)

p <- plot_mnsl(col)
summary(p)

p + ggplot2::facet_wrap(~ num, nrow = 2)

p+ggplot2::facet_wrap(~soil_col$munsell[4:46], nc = 6)
?plot_mnsl

sample()

#soil
#git