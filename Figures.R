#The figures

#Grid ####

x <- rep(seq(0,50, 10), 6)
y <- c(rep(0,6), rep(10,6), rep(20,6), rep(30,6), rep(40, 6), rep(50,6))
x2 <- c(5, 15, 0, 25, 30,50,35,5,45,45,15,10,40,15)
y2 <- c(0, 45, 35, 30, 35, 5, 40, 20, 15, 50, 10, 45, 25, 15)


plot(x,y,
     pch = 0, 
     xaxt = "n", 
     yaxt = "n", 
     xlab = "", 
     ylab = "", 
     frame.plot = F, 
     cex = 2.3,
     lwd = 2,
     panel.first = abline(25, 0, lwd = 950, col = "gray93"))

axis(side = 1, line = 1, at = c(0,50), labels = c("", ""))
axis(side = 2, line = 0.5, at = c(20,30), labels = c("", ""))
title(ylab = "10 m ", line = 1.8, cex.lab = 1.2)
title(xlab = "50 m ", line = 1.5, cex.lab = 1.2)



points(x2, y2,
       pch = 0,
       cex = 2.3,
       col = "blue",
       lwd = 2)




#Figure 2 ####

load("occ_model_scaled_nospat_17Aug.rda")
var <- calc.varpart(occ_model_scaled_nospat_17Aug)
part <- data.frame(enviro = var$varpart.X, bio = var$varpart.lv)
bardat <- data.frame(Species =rep(rownames(part),2), variance = c(as.vector(part[,1]), as.vector(part[,2])), type = c(rep("enviro",10), rep("lv",10)))

(occ <- ggplot(bardat, aes(x=Species, y=variance, fill = type)) +
  geom_bar(stat="identity", width = 0.55, position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 40, 
                                   colour = "black", 
                                   size = 13,
                                   vjust = 0.98, 
                                   hjust = 1,
                                   face = "italic"),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11,
                                   colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  ylab("Proportion of variance") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.001)) +
  scale_x_discrete(
    limits = c("R_comptonii",
               "R_burtoniae", 
               "Oophytum_sp",
               "Dicrocaulon_sp",
               "D_diversifolium",
               "C_staminodiosum",
               "A_delaetii",
               "A_framesii",
               "C_spissum",
               "A_fissum"),
    
    labels = c("R. comptonii",
               "R. burtoniae", 
               "Oophytum sp.",
               "Dicrocaulon sp.",
               "D. diversifolium",
               "C. staminodiosum",
               "A. delaetii",
               "A. framesii",
               "C. spissum",
               "A. fissum"
  )) +
  scale_fill_manual(labels = c("Environment", "Latent"),
                    name = "Source",
                    values = c("#F26419", "#F6B540")) +
  guides(fill = F) +
  annotate("text", x=1, y=0.99, label= "b)", size = 10))

#abundance

load("poisson_abn_1Sept.rda")
var2 <- calc.varpart(poisson_abn_1Sept)
part <- data.frame(enviro = var2$varpart.X, bio = var2$varpart.lv)
bardat2 <- data.frame(Species =rep(rownames(part),2), variance = c(as.vector(part[,1]), as.vector(part[,2])), type = c(rep("enviro",10), rep("lv",10)))

(abn <- ggplot(bardat2, aes(x=Species, y=variance, fill = type)) +
  geom_bar(stat="identity", width = 0.55, position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   colour = "black", 
                                   size = 10,
                                   vjust = 0.98, 
                                   hjust = 1,
                                   face = "italic"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 11,
                                   colour = "black"),
        panel.grid = element_blank(),
        legend.position=c(0.85, 0.9),
        legend.text = element_text(size = 15)) +
  ylab("Proportion of variance") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.001)) +
  scale_x_discrete(
    limits = c("R_comptonii",
               "R_burtoniae", 
               "Oophytum_sp",
               "Dicrocaulon_sp",
               "D_diversifolium",
               "C_staminodiosum",
               "A_delaetii",
               "A_framesii",
               "C_spissum",
               "A_fissum"),
    
    labels = NULL) + scale_fill_manual(labels = c("Environment", "Latent"),
                          name = "",
                          values = c("#F26419", "#F6B540")
                    ) +
  annotate("text", x=1, y=0.99, label= "a)", size = 10))
  

cowplot::plot_grid(abn,occ, ncol = 1, rel_heights = c(0.82,1))
  
  
  
  #Figure 3 ####
load("occ_model_scaled_nospat_17Aug.rda")

envcors <- get.enviro.cor(occ_model_scaled_nospat_17Aug)
rescors <- get.residual.cor(occ_model_scaled_nospat_17Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor, order = "hclust")
corrplot(rescors$sig.correlaton, order = "AOE", type = "lower")
corrplot(rescors$precision)

par(mfrow = c(1,2))
qgraph::qgraph(envcors$sig.cor, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10, labels = c("Rbr", "Rcp", "Ddv", "Adl", "Afr", "Afs", "Csp", "Cst", "Dcr", "Oph"))
text(-1,1, "a)", cex = 1.5)
qgraph::qgraph(rescors$sig.correlaton, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10
, labels = c("Rbr", "Rcp", "Ddv", "Adl", "Afr", "Afs", "Csp", "Cst", "Dcr", "Oph"))
text(-1,1, "b)", cex = 1.5)

env <- as.dist(envcors$cor)
res <- as.dist(rescors$correlation) 
plot(env, res)
cor.test(env, res)


plot(envcors$cor, rescors$correlation,
     xlim = c(-1,1))
abline(0,0)
abline(v = 0)

mod <- lm(rescors$correlation ~ envcors$cor)
abline(mod)
summary(mod)
cor.test(rescors$correlation, envcors$cor)

require("lessR")

cormat <- as.data.frame(envcors$sig.cor)
cormat <- as.data.frame(rescors$correlation)

names(cormat) <- c("R. burtoniae", 
                   "R. comptonii", 
                   "D. diversifolium",
                   "A. delaetii",
                   "A. fissum",
                   "A. framesii",
                   "C. spissum",
                   "C. staminodiosum",
                   "Dicrocaulon sp.",
                   "Oophytum sp.")
row.names(cormat) <- c("R. burtoniae", 
                                        "R. comptonii", 
                                        "D. diversifolium",
                                        "A. delaetii",
                                        "A. fissum",
                                        "A. framesii",
                                        "C. spissum",
                                        "C. staminodiosum",
                                        "Dicrocaulon sp.",
                                        "Oophytum sp.")



cormat <- cormat %>% select("Oophytum sp.",
                   "A. delaetii",
                   "A. framesii",
                   "C. staminodiosum",
                   "C. spissum",
                   "A. fissum",
                   "R. comptonii",
                   "R. burtoniae", 
                   "D. diversifolium",
                   "Dicrocaulon sp.")

new <- cormat[ colnames(cormat), ] %>% as.matrix()

corrplot(new,
         #col = brewer.pal(n = 8, name = "RdYlBu"),
        
         tl.col = "black",
         tl.cex = 1.1,
         font = 3,
         type = "upper")

text


qgraph::qgraph(new, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10,
               labels = c( "Oph", "Adl", "Afr", "Cst","Csp","Afs", "Rcp", "Rbr","Ddv","Dcr" ))


#Figure 5 ####
load("poisson_abn_1Sept.rda")
mod <- occ_model_scaled_nospat_17Aug
mod <- poisson_abn_1Sept

xlabs <-  c(
  "R. burtoniae", 
  "R. comptonii", 
  "D. diversifolium",
  "A. delaetii",
  "A. fissum",
  "A. framesii",
  "C. spissum",
  "C. staminodiosum",
  "Dicrocaulon sp.",
  "Oophytum sp.")

ylabs <- c("pH",
           "Na",
           "C:N ratio",
           "Ca",
           "quartz",
           "elevation",
           "%1-2mm",
           "%2-5mm",
           "P",
           "%N",
           "aspect",
           "clay",
           expression(paste(delta^15~N)),
           "drainage",
           "slope")

par(mfrow = c(1,3))


#R burt ####
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[1, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[1, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[1, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[1, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[1],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
  
)

axis(
  2,
  labels = ylabs,
  at = 1:15,
  las = 2,
  cex.axis = 1.7
)

segments(
  x0 = mod$hpdintervals$X.coefs[1, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[1, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#R comp ####

par(mar = c(5.1, 4.5, 2.1, 2.5) )
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[2, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[2, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[2, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[2, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[2],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[2, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[2, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#D div ####
#par(mar = c(5.1, 3, 2.1, 2))
#par(mar = c(5.1, 6.5, 2.1, 0))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[3, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[3, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[3, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[3, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[3],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[3, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[3, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#A del ####

#par(mar = c(5.1, 2, 2.1, 4.5))
par(mar = c(5.1, 7, 2.1, 0))
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[4, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[4, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[4, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[4, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[4],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[4, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[4, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

axis(
  2,
  labels = ylabs,
  at = 1:15,
  las = 2,
  cex.axis = 1.7
)


#A fiss ####
#par(mar = c(5.1, 1, 2.1, 4))
par(mar = c(5.1, 4.5, 2.1, 2.5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[5, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[5, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[5, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[5, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[5],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[5, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[5, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fra ####

#par(mar = c(5.1, 5, 2.1, 0))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[6, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[6, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[6, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[6, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[6],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[6, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[6, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#C spis ####
#par(mar = c(5.1, 4, 2.1, 1) )
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[7, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[7, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[7, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[7, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[7],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)
axis(
  2,
  labels = ylabs,
  at = 1:15,
  las = 2,
  cex.axis = 1.7
)


segments(
  x0 = mod$hpdintervals$X.coefs[7, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[7, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#C sta ####

#par(mar = c(5.1, 3, 2.1, 2))
par(mar = c(5.1, 4.5, 2.1, 2.5) )
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[8, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[8, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[8, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[8, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[8],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[8, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[8, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#Dicr ####
#par(mar = c(5.1, 2, 2.1, 3))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[9, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[9, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[9, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[9, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[9],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[9, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[9, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#Ooph ####

#par(mar = c(5.1, 1, 2.1, 4))
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[10, 1:15, "lower"]))
col.seq[mod$hpdintervals$X.coefs[10, 1:15, "lower"] < 0 &
          mod$hpdintervals$X.coefs[10, 1:15, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[10, 1:15]),
  y = 1:15,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[10],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:15, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:15, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[10, 1:15, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[10, 1:15, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

axis(
  2,
  labels = ylabs,
  at = 1:15,
  las = 2,
  cex.axis = 1.7
)

#Occurrence est ####
load("occ_model_scaled_nospat_17Aug.rda")
mod <- occ_model_scaled_nospat_17Aug


xlabs <-  c(
  "R. burtoniae", 
  "R. comptonii", 
  "D. diversifolium",
  "A. delaetii",
  "A. fissum",
  "A. framesii",
  "C. spissum",
  "C. staminodiosum",
  "Dicrocaulon sp.",
  "Oophytum sp.")

ylabs <- c("pH",
           "Na",
           "C:N ratio",
           "Ca",
           "quartz",
           "elevation",
           "%1-2mm",
           "%2-5mm",
           "P",
           "%N",
           "aspect",
           "clay",
           expression(paste(delta^15~N)),
           "drainage")

par(mfrow = c(1,3))


#R burt ####
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[1, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[1, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[1, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[1, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[1],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
  
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.7
)

segments(
  x0 = mod$hpdintervals$X.coefs[1, 1:14, "lower"],
  y0 = 1:15,
  x1 = mod$hpdintervals$X.coefs[1, 1:14, "upper"],
  y1 = 1:15,
  col = col.seq
)
abline(v = 0, lty = 3)

#R comp ####

par(mar = c(5.1, 4.5, 2.1, 2.5) )
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[2, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[2, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[2, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[2, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[2],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[2, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[2, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#D div ####
#par(mar = c(5.1, 3, 2.1, 2))
#par(mar = c(5.1, 6.5, 2.1, 0))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[3, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[3, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[3, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[3, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[3],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[3, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[3, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A del ####

#par(mar = c(5.1, 2, 2.1, 4.5))
par(mar = c(5.1, 7, 2.1, 0))
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[4, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[4, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[4, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[4, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[4],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[4, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[4, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.7
)


#A fiss ####
#par(mar = c(5.1, 1, 2.1, 4))
par(mar = c(5.1, 4.5, 2.1, 2.5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[5, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[5, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[5, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[5, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[5],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[5, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[5, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fra ####

#par(mar = c(5.1, 5, 2.1, 0))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[6, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[6, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[6, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[6, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[6],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[6, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[6, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#C spis ####
#par(mar = c(5.1, 4, 2.1, 1) )
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[7, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[7, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[7, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[7, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[7],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)
axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.7
)


segments(
  x0 = mod$hpdintervals$X.coefs[7, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[7, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#C sta ####

#par(mar = c(5.1, 3, 2.1, 2))
par(mar = c(5.1, 4.5, 2.1, 2.5) )
col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[8, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[8, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[8, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[8, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[8],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[8, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[8, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#Dicr ####
#par(mar = c(5.1, 2, 2.1, 3))
par(mar = c(5.1, 2, 2.1, 5) )

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[9, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[9, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[9, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[9, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[9],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)

#axis(
#  2,
#  labels = ylabs,
#  at = 1:15,
#  las = 2,
#  cex.axis = 1.5
#)

segments(
  x0 = mod$hpdintervals$X.coefs[9, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[9, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#Ooph ####

#par(mar = c(5.1, 1, 2.1, 4))
par(mar = c(5.1, 7, 2.1, 0))

col.seq <-
  rep("black", length(mod$hpdintervals$X.coefs[10, 1:14, "lower"]))
col.seq[mod$hpdintervals$X.coefs[10, 1:14, "lower"] < 0 &
          mod$hpdintervals$X.coefs[10, 1:14, "upper"] > 0] <- "grey"

plot(
  x = c(mod$X.coefs.mean[10, 1:14]),
  y = 1:14,
  yaxt = "n",
  ylab = "",
  xlab = xlabs[10],
  font.lab = 3,
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex = 1.5,
  cex.axis = 1.6,
  cex.lab = 2.3
)


segments(
  x0 = mod$hpdintervals$X.coefs[10, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[10, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.7
)






#Spatial ####

load("/Users/larawootton/Documents/Honours/Data/plots12_19Aug.rda")
test <- all_dat %>% filter(site == "site3")
newpred1 <- predict.boral(plots12_19Aug, 
                          newX = test[,13:26], 
                          predict.type = "marginal",
                          est = "mean")
pred <- newpred1$linpred
pred_all1 <- ROCR::prediction(pred[,c(1:5,7:10)], test[,c(1:5,7:10)])
aucs1 <- ROCR::performance(pred_all1, "auc")
au <- unlist(aucs1@y.values)
Mod1 <- c(au[1:5], NA, au[6:9])

#Model 2

load("/Users/larawootton/Documents/Honours/Data/plots13_1Aug.rda")
test <- all_dat %>% filter(site == "site2")
newpred2 <- predict.boral(plots13_19Aug, 
                          newX = test[,13:26], 
                          predict.type = "marginal",
                          est = "mean")

pred <- newpred2$linpred
pred_all2 <- ROCR::prediction(pred[,c(1, 3:7,9:10)], test[,c(1, 3:7,9:10)])
aucs2 <- ROCR::performance(pred_all2, "auc")
au <- unlist(aucs2@y.values)
Mod2 <- c(au[1], NA, au[2:6], NA, au[7:8])

#Model 3

load("/Users/larawootton/Documents/Honours/Data/plots23_19Aug.rda")
test <- all_dat %>% filter(site == "site1")
newpred3 <- predict.boral(plots23_19Aug, 
                          newX = test[,13:26], 
                          predict.type = "marginal",
                          est = "mean")

pred <- newpred3$linpred
pred_all3 <- ROCR::prediction(pred, test[,1:10])
aucs3 <- ROCR::performance(pred_all3, "auc")
Mod3 <- unlist(aucs3@y.values)

#Model 4

load("grid_rand_occ_29Aug_scale.rda")
typedat <- cbind(all_dat, type = type)
grid <- typedat %>% filter(type == "grid") %>% select(-type)
test <- typedat %>% filter(type == "random") %>% select(-type)

newpred6 <- predict.boral(grid_rand_occ_29Aug_scale, 
                          newX = test[,13:26], 
                          predict.type = "marginal",
                          est = "mean")
pred <- newpred6$linpred
pred_all6 <- ROCR::prediction(pred, test[,1:10])
aucs6 <- ROCR::performance(pred_all6, "auc")
Mod4 <- unlist(aucs6@y.values)

mods <- data.frame(auc = c(Mod1, Mod2, Mod3, Mod4), 
                   model = c(rep("Model 1", 10),rep("Model 2", 10),rep("Model 3", 10),rep("Model 4", 10)))

boxplot(mods$auc ~ mods$model,
        whisklty = 1,
        staplelty = 0,
        boxwex = 0.3,
        ylab = "AUC")
points(c(mean(Mod1, na.rm = T), mean(Mod2, na.rm = T), mean(Mod3), mean(Mod4)), 
         pch = 19,
         col = "red",
       cex = 1.4)
mtext(c("ab", "a", "ab", "b"), at = c(1:4), cex = 1.5)


ano <- aov(mods$auc ~ mods$model)
summary(ano) #F = 3.359, p = 0.0303
tuk <- TukeyHSD(ano)
plot(tuk)
summary(tuk)

#By species

sp_mods <- data.frame(species = colnames(all_dat[,1:10]), Mod1, Mod2, Mod3)



boxplot(mods[1,1])


#Visualising the models
my_dat <- my_soil_df %>% filter(!is.na(type)) %>% na.omit()

par(mfrow = c(2,2))
plot(my_dat$lon, my_dat$lat, 
     pch = 22, cex = 1.4, 
     col = ifelse(my_dat$site == "site1" | my_dat$site=="site2", "red", "blue"),
     main = "Model 1: Sites 1 and 2 predicting to site 3",
     ylab = "Lat",
     xlab = "",
     cex.lab = 1.4)

plot(my_dat$lon, my_dat$lat, 
     pch = 22, cex = 1.4, 
     col = ifelse(my_dat$site == "site3" | my_dat$site=="site1", "red", "blue"),
     main = "Model 2: Sites 1 and 3 predicting to site 2",
     ylab = "",
     xlab = "",
     cex.lab = 1.4)

plot(my_dat$lon, my_dat$lat, 
     pch = 22, cex = 1.4, 
     col = ifelse(my_dat$site == "site3" | my_dat$site=="site2", "red", "blue"),
     main = "Model 3: Sites 2 and 3 predicting to site 1",
     ylab = "Lat",
     xlab = "Lon",
     cex.lab = 1.4)

plot(my_dat$lon, my_dat$lat, 
     pch = 22, cex = 1.4, 
     col = ifelse(my_dat$type == "grid", "red", "blue"),
     main = "Model 4: Use grid plots to predict onto random plots",
     ylab = "",
     xlab = "Lon",
     cex.lab = 1.4)



##Tables####

#Soil variables####

N_percent <- my_soil_df %>% 
  filter(!is.na(type)) %>% 
  summarise(Average = mean(N_perc, na.rm = T),
                                      SD = sd(N_perc, na.rm = T),
                                      Median = median(N_perc, na.rm = T),
                                      Max = max(N_perc, na.rm = T),
                                      Min = min(N_perc, na.rm = T))
 

dN_corrected <-  my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(corr_dN, na.rm = T),
                                          SD = sd(corr_dN, na.rm = T),
                                          Median = median(corr_dN, na.rm = T),
                                          Max = max(corr_dN, na.rm = T),
                                          Min = min(corr_dN, na.rm = T))

C_percent <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(C_perc, na.rm = T),
                                      SD = sd(C_perc, na.rm = T),
                                      Median = median(C_perc, na.rm = T),
                                      Max = max(C_perc, na.rm = T),
                                      Min = min(C_perc, na.rm = T))
dC <-  my_soil_df %>% summarise(Average = mean(corr_dC, na.rm = T),
                                SD = sd(corr_dC, na.rm = T),
                                Median = median(corr_dC, na.rm = T),
                                Max = max(corr_dC, na.rm = T),
                                Min = min(dC, na.rm = T))


C_N_ratio <-  my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(C_N_ratio, na.rm = T),
                                       SD = sd(C_N_ratio, na.rm = T),
                                       Median = median(C_N_ratio, na.rm = T),
                                       Max = max(C_N_ratio, na.rm = T),
                                       Min = min(C_N_ratio, na.rm = T))
Acidity <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(acidity),
                                  SD = sd(acidity),
                                  Median = median(acidity),
                                  Max = max(acidity),
                                  Min = min(acidity))
Ca <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(Ca),
                             SD = sd(Ca),
                             Median = median(Ca),
                             Max = max(Ca),
                             Min = min(Ca))
Mg <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(Mg),
                             SD = sd(Mg),
                             Median = median(Mg),
                             Max = max(Mg),
                             Min = min(Mg))
Na <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(Na),
                             SD = sd(Na),
                             Median = median(Na),
                             Max = max(Na),
                             Min = min(Na))
K <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(K),
                            SD = sd(K),
                            Median = median(P),
                            Max = max(K),
                            Min = min(K))
P <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(P),
                            SD = sd(P),
                            Median = median(P),
                            Max = max(P),
                            Min = min(P))
Olsen <- my_soil_df %>% 
  filter(!is.na(type)) %>%
  summarise(Average = mean(Olsen),
                                SD = sd(Olsen),
                                Median = median(Olsen),
                                Max = max(Olsen),
                                Min = min(Olsen))
clay <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(Clay, na.rm = T), 
                                     SD = sd(Clay, na.rm = T),
                                     Median = median(Clay, na.rm = T),
                              Max = max(Clay, na.rm = T), 
                              Min = min(Clay, na.rm = T) 
                              )
silt <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(Silt, na.rm = T),
                                     SD = sd(Silt, na.rm = T),
                                     Median = median(Silt, na.rm = T),
                              Max = max(Silt, na.rm = T), 
                              Min = min(Silt, na.rm = T)
                              ) 
total_sand <-my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(Sand, na.rm = T), 
                                     SD = sd(Sand, na.rm = T),
                                     Median = median(Sand, na.rm = T),
                                    Max = max(Sand, na.rm = T), 
                                    Min = min(Sand, na.rm = T) 
                                    )
less2 <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(percent_over1, na.rm = T),
                                     SD = sd(percent_over1, na.rm = T),
                                     Median = median(percent_over1, na.rm = T),
                               Max = max(percent_over1, na.rm = T), 
                               Min = min(percent_over1, na.rm = T) 
                               )
less5 <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(percent_over2, na.rm = T),
                                     SD = sd(percent_over2, na.rm = T),
                                     Median = median(percent_over2, na.rm = T),
                               Max = max(percent_over2, na.rm = T), 
                               Min = min(percent_over2, na.rm = T)
                              )

cond <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(conductivity_ms, na.rm = T),
                                     
                              SD = sd(conductivity_ms, na.rm = T),
                              Median = median(conductivity_ms),
                              Min = min(conductivity_ms, na.rm = T),
                              Max = max(conductivity_ms, na.rm = T))
ph_kcl <- my_soil_df %>% 
  filter(!is.na(type)) %>% summarise(Average = mean(ph_kcl, na.rm = T),
                                    
                                SD = sd(ph_kcl, na.rm = T),
                                Median = median(ph_kcl),
                                Min = min(ph_kcl, na.rm = T),
                                Max = max(ph_kcl, na.rm = T))

measure <- rbind(Acidity,
      ph_kcl,
      cond,
      Na,
      Ca, 
      Mg,
      K,
      P,
      Olsen,
      dC,
      C_percent,
      dN_corrected,
      N_percent,
      C_N_ratio,
      clay,
      silt,
      total_sand,
      less2,
      less5) %>% round(2)

write.csv(measure, file = "Soil_table.csv", row.names = F)


