#!/usr/bin/env Rscript

source("scripts/plotpredict.R")
source("scripts/summarypop.R")


source("scripts/Gmatrices.R")


tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")



pdf("figure2.pdf", width=12, height=8)
layout(cbind(1:2,3:4))
par(mar=c(3, 4.2, 3, 1))
plot.updown(tovar, tovar.mat$G, tovar.mat$Gv, tovar.mat$P, main="Tovar", ylab=expression("log GA (mm"^2*")"), xlab="", ylim=c(-0.4,0.4))
par(mar=c(5, 4.2, 1, 1))
plot.updown(tovar, tovar.mat$G, tovar.mat$Gv, tovar.mat$P, target="mu.y", main="", ylab=expression("log GBA (mm"^2*")"), ylim=c(-0.3,0.3))
par(mar=c(3, 4.2, 3, 1))
plot.updown(tulum, tulum.mat$G, tulum.mat$Gv, tulum.mat$P, main="Tulum", ylab="", xlab="", ylim=c(-0.4,0.4))
par(mar=c(5, 4.2, 1, 1))
plot.updown(tulum, tulum.mat$G, tulum.mat$Gv, tulum.mat$P, target="mu.y", main="", ylab="", ylim=c(-0.3,0.3))
dev.off()
