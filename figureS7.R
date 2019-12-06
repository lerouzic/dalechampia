#!/usr/bin/env Rscript

source("scripts/plotpredict.R")


source("data/matrices.R") # G, P, and Gv matrices
tovar <- read.table("data/tovar.txt", header=TRUE, na.strings="-")
tulum <- read.table("data/tulum.txt", header=TRUE, na.strings="-")



pdf("figureS7.pdf", width=12, height=8)
layout(cbind(1:2,3:4))
par(mar=c(3, 4.2, 3, 1))
plot.raw(tovar, G.tovar, Gv.tovar, P.tovar, main="Tovar", ylab=expression("log GA (mm"^2*")"), xlab="", ylim=c(2.4,3.7))
par(mar=c(5, 4.2, 1, 1))
plot.raw(tovar, G.tovar, Gv.tovar, P.tovar, target="mu.y", main="", ylab=expression("log GBA (mm"^2*")"), ylim=c(5.6,6.4))
par(mar=c(3, 4.2, 3, 1))
plot.raw(tulum, G.tulum, Gv.tulum, P.tulum, main="Tulum", ylab="", xlab="", ylim=c(2.4,3.7))
par(mar=c(5, 4.2, 1, 1))
plot.raw(tulum, G.tulum, Gv.tulum, P.tulum, target="mu.y", main="", ylab="", ylim=c(5.6,6.4))
dev.off()
