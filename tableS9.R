#!/usr/bin/env Rscript

source("scripts/plotpredict.R")


source("data/matrices.R") # G, P, and Gv matrices
tovar <- read.table("data/tovar.txt", header=TRUE, na.strings="-")
tulum <- read.table("data/tulum.txt", header=TRUE, na.strings="-")

# Tovar, GA
toga <- cbind(
	plot.raw(tovar, G.tovar, Gv.tovar, P.tovar, verr.only=TRUE), 
	vControl=plot.control(tovar, G.tovar, Gv.tovar, P.tovar, verr.only=TRUE)$Vtot, 
	Vselect=plot.updown(tovar, G.tovar, Gv.tovar, P.tovar, verr.only=TRUE)$Vtot)
write.table(toga, "tableS8-TovGA.txt", quote=FALSE, sep="\t", row.names=FALSE)


# Tovar, UBA
touba <- cbind(
	plot.raw(tovar, G.tovar, Gv.tovar, P.tovar, target="mu.y", verr.only=TRUE), 
	vControl=plot.control(tovar, G.tovar, Gv.tovar, P.tovar, target="mu.y", verr.only=TRUE)$Vtot, 
	Vselect=plot.updown(tovar, G.tovar, Gv.tovar, P.tovar, target="mu.y", verr.only=TRUE)$Vtot)
write.table(touba, "tableS8-TovUBA.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Tulum, GA
tuga <- cbind(
	plot.raw(tulum, G.tulum, Gv.tulum, P.tulum, verr.only=TRUE), 
	vControl=plot.control(tulum, G.tulum, Gv.tulum, P.tulum, verr.only=TRUE)$Vtot, 
	Vselect=plot.updown(tulum, G.tulum, Gv.tulum, P.tulum, verr.only=TRUE)$Vtot)
write.table(tuga, "tableS8-TulGA.txt", quote=FALSE, sep="\t", row.names=FALSE)


# Tulum, UBA
tuuba <- cbind(
	plot.raw(tulum, G.tulum, Gv.tulum, P.tulum, target="mu.y", verr.only=TRUE), 
	vControl=plot.control(tulum, G.tulum, Gv.tulum, P.tulum, target="mu.y", verr.only=TRUE)$Vtot, 
	Vselect=plot.updown(tulum, G.tulum, Gv.tulum, P.tulum, target="mu.y", verr.only=TRUE)$Vtot)
write.table(tuuba, "tableS8-tulUBA.txt", quote=FALSE, sep="\t", row.names=FALSE)
