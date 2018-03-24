
#Creating a figure showing the experimental set up  

x1 <- as.vector(c(rep(1,11), rep(2,11), rep(3,11), rep(4,11), rep(5,11), rep(6,11), rep(7,11), rep(8,11), rep(9,11), rep(10,11), rep(11,11)))
y1 <- as.vector(rep(seq(1, 11, 1), 11))


x3 <- as.vector(c(rep(1,6), rep(3,6), rep(5,6), rep(7,6), rep(9,6), rep(11,6)))
y3 <- as.vector(rep(seq(1, 11, 2), 6))

dat <- as.data.frame(cbind(x1, y1))
dat2 <- data.frame(x3,y3)

plot(dat$x1, dat$y1)
plot(dat2$x2, dat2$y2)

plot(dat$x1, dat$y1, 
     pch = 19, xaxt = "n", yaxt = "n", xlab = "", ylab = "", frame.plot = F, cex = 1.2)
points(dat2$x3, dat2$y3, cex = 4.2)

axis(side = 1, line = 1, at = c(1,11), labels = c("", ""))
axis(side = 2, line = 1, at = c(1,11), labels = c("", ""))
title(ylab = "50 m ", line = 1.2, cex.lab = 1.2)
title(xlab = "50 m ", line = 1.2, cex.lab = 1.2)

axis(side = 3, line = 1, labels = c("0", "5", "10"))