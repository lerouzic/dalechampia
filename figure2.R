#!/usr/bin/env Rscript

source("scripts/plotpredict.R")
source("scripts/summarypop.R")


source("scripts/Gmatrices.R")


tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")



pdf("figure2.pdf", width=8, height=8)
layout(cbind(1:2,3:4))
par(mar=c(3, 4.2, 3, 1))
plot.data.recenter(
	recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.x", normalization="updown"),
	main="Tovar", ylab=expression("log GA (mm"^2*")"), xlab="", ylim=c(-0.4,0.4))
par(mar=c(5, 4.2, 1, 1))
plot.data.recenter(
	recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.y", normalization="updown"),
	main="", ylab=expression("log UBA (mm"^2*")"), xlab="Generation", ylim=c(-0.3,0.3))
par(mar=c(3, 4.2, 3, 1))
plot.data.recenter(
	recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.x", normalization="updown"),
	main="Tulum", ylab="", xlab="", ylim=c(-0.4,0.4))
par(mar=c(5, 4.2, 1, 1))
plot.data.recenter(
	recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.y", normalization="updown"),
	main="", ylab="", xlab="Generation", ylim=c(-0.3,0.3))
dev.off()
