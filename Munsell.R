require("munsell")

col <- c("10YR 4/6",
"10YR 5/6",
"10YR 6/6",
"10YR 5/8",
"10YR 6/6",
"10YR 6/6",
"10YR 5/4",
"10YR 4/6",
"7.5YR 5/6",
"10YR 6/6",
"10YR 6/6",
"10YR 6/6",
"7.5YR 5/6",
"10YR 5/4",
"10YR 6/4",
"10YR 6/6",
"10YR 6/4")

setwd("/Users/larawootton/Documents/Honours")
soil_col <- read.csv("Munsell_8Mar.csv", sep = ";")
grid <- subset(soil_col, type == "grid")


col <- soil_col$munsell[4:39]

plot_mnsl(col)
plot_mnsl(grid$munsell[1:20], na.rm = T)

p <- plot_mnsl
summary(p)

p+ggplot2::facet_wrap(~soil_col$munsell[4:46], ncol = 6)
?plot_mnsl

sample()

#soil
#git