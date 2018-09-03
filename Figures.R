#The figures


load("occ_model_scaled_nospat_17Aug.rda")

#Figure 2 ####
var <- calc.varpart(occ_model_scaled_nospat_17Aug)
part <- data.frame(enviro = var$varpart.X, bio = var$varpart.lv)
bardat <- data.frame(Species =rep(rownames(part),2), variance = c(as.vector(part[,1]), as.vector(part[,2])), type = c(rep("enviro",10), rep("lv",10)))

ggplot(bardat, aes(x=Species, y=variance, fill = type)) +
  geom_bar(stat="identity", width = 0.55, position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   colour = "black", 
                                   size = 10,
                                   vjust = 0.98, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  ylab("Proportion of variance") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.001)) +
  scale_x_discrete(
    limits = c(
      "Oophytum_sp",
      "C_staminodiosum",
      "R_burtoniae",
      "A_framesii",
      "D_diversifolium",
      "R_comptonii",
      "Dicrocaulon_sp",
      "A_delaetii",
      "A_fissum",
      "C_spissum"),
    labels = c("Oophytum sp.",
               "C. staminodiosum",
               "R. burtoniae", 
               "A. framesii",
               "D. diversifolium",
               "R. comptonii",
               "Dicrocaulon sp.",
               "A. delaetii",
               "A. fissum",
               "C. spissum"
    )
  ) +
  scale_fill_manual(labels = c("Environment", "Latent"),
                    name = "Source",
                    values = c("#F26419", "#F6B540"))

#Figure 3 ####

envcors <- get.enviro.cor(occ_model_scaled_nospat_17Aug)
rescors <- get.residual.cor(occ_model_scaled_nospat_17Aug)
corrplot(envcors$cor)
corrplot(envcors$sig.cor, order = "hclust")
corrplot(rescors$correlation, order = "AOE", type = "lower")
corrplot(rescors$sig.correlaton)

par(mfrow = c(1,2))
qgraph::qgraph(envcors$sig.cor, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10, labels = c("Rbr", "Rcp", "Ddv", "Adl", "Afr", "Afs", "Cst", "Csp", "Dcr", "Oph"))
text(-1,1, "a)", cex = 1.5)
qgraph::qgraph(rescors$sig.correlaton, shape="circle", posCol="darkgreen", negCol="darkred", layout="groups", vsize=10, labels = c("Rbr", "Rcp", "Ddv", "Adl", "Afr", "Afs", "Cst", "Csp", "Dcr", "Oph"))
text(-1,1, "b)", cex = 1.5)

#Figure 4 ####

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

ylabs <- c("pH et al",
           "salt",
           "carbon",
           "Ca",
           "quartz",
           "elevation",
           "%1-2mm",
           "%2-5mm",
           "P",
           "%N",
           "aspect",
           "texture",
           "d15N",
           "drainage")
par(mfrow = c(5,2))


#R burt
par(mar = c(2.1, 5, 2.1, 1.1))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1
)

segments(
  x0 = mod$hpdintervals$X.coefs[1, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[1, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#R comp

par(mar = c(2.1, 1.1, 2.1, 5))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)


segments(
  x0 = mod$hpdintervals$X.coefs[2, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[2, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#D div
par(mar = c(2.1, 5, 2.1, 1.1))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1
)

segments(
  x0 = mod$hpdintervals$X.coefs[3, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[3, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A del

par(mar = c(2.1, 1.1, 2.1, 5))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)


segments(
  x0 = mod$hpdintervals$X.coefs[4, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[4, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fiss
par(mar = c(2.1, 5, 2.1, 1.1))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1
)

segments(
  x0 = mod$hpdintervals$X.coefs[5, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[5, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fra

par(mar = c(2.1, 1.1, 2.1, 5))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)


segments(
  x0 = mod$hpdintervals$X.coefs[6, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[6, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#C sta
par(mar = c(2.1, 5, 2.1, 1.1))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1
)

segments(
  x0 = mod$hpdintervals$X.coefs[7, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[7, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#C spis

par(mar = c(2.1, 1.1, 2.1, 5))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)


segments(
  x0 = mod$hpdintervals$X.coefs[8, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[8, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#Dicr
par(mar = c(2.1, 5, 2.1, 1.1))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1
)

segments(
  x0 = mod$hpdintervals$X.coefs[9, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[9, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#Ooph

par(mar = c(4.1, 1.1, 2.1, 5))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x"
)


segments(
  x0 = mod$hpdintervals$X.coefs[10, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[10, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)


#2x5####

par(mfrow = c(2,5))


#R burt
par(mar = c(5.1, 5, 2.1, 0))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.1
)

segments(
  x0 = mod$hpdintervals$X.coefs[1, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[1, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#R comp ####

par(mar = c(5.1, 4, 2.1, 1) )
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
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
par(mar = c(5.1, 3, 2.1, 2))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)


segments(
  x0 = mod$hpdintervals$X.coefs[3, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[3, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A del ####

par(mar = c(5.1, 2, 2.1, 3))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)


segments(
  x0 = mod$hpdintervals$X.coefs[4, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[4, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fiss ####
par(mar = c(5.1, 1, 2.1, 4))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)



segments(
  x0 = mod$hpdintervals$X.coefs[5, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[5, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#A fra ####

par(mar = c(5.1, 5, 2.1, 0))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)

axis(
  2,
  labels = ylabs,
  at = 1:14,
  las = 2,
  cex.axis = 1.1
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
par(mar = c(5.1, 4, 2.1, 1) )

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
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

par(mar = c(5.1, 3, 2.1, 2))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
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
par(mar = c(5.1, 2, 2.1, 3))

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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)


segments(
  x0 = mod$hpdintervals$X.coefs[9, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[9, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)

#Ooph

par(mar = c(5.1, 1, 2.1, 4))
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
  col = col.seq,
  xlim = c(
    min(mod$hpdintervals$X.coefs[, 1:14, "lower"]),
    max(mod$hpdintervals$X.coefs[, 1:14, "upper"])
  ),
  pch = "x",
  cex.lab = 1.4
)


segments(
  x0 = mod$hpdintervals$X.coefs[10, 1:14, "lower"],
  y0 = 1:14,
  x1 = mod$hpdintervals$X.coefs[10, 1:14, "upper"],
  y1 = 1:14,
  col = col.seq
)
abline(v = 0, lty = 3)
