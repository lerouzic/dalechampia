#!/usr/bin/env Rscript

# Table 1. G matrices are estimated from an animal model fit on the diallel data

source("scripts/Gmatrices.R", chdir=TRUE)

dd <- 2 # number of digits for rounding
mf <- 100 # multiplicaton factor

for (excluding.self in c(FALSE, TRUE)) {
	tovar <- get.matrices("Tovar", excluding.self=excluding.self)
	tulum <- get.matrices("Tulum", excluding.self=excluding.self)
	
	tovar.G.cor <- tovar$G.VCV[,2]/sqrt(tovar$G.VCV[,1])/sqrt(tovar$G.VCV[,4])
	tovar.E.cor <- tovar$E.VCV[,2]/sqrt(tovar$E.VCV[,1])/sqrt(tovar$E.VCV[,4])
	tovar.P.cor <- tovar$P.VCV[,2]/sqrt(tovar$P.VCV[,1])/sqrt(tovar$P.VCV[,4])
	
	tulum.G.cor <- tulum$G.VCV[,2]/sqrt(tulum$G.VCV[,1])/sqrt(tulum$G.VCV[,4])
	tulum.E.cor <- tulum$E.VCV[,2]/sqrt(tulum$E.VCV[,1])/sqrt(tulum$E.VCV[,4])
	tulum.P.cor <- tulum$P.VCV[,2]/sqrt(tulum$P.VCV[,1])/sqrt(tulum$P.VCV[,4])
	
	
	data.tab1.col1 <- list(tovar$G.VCV[,1], tovar.G.cor, tovar$E.VCV[,1], tovar.E.cor, tovar$P.VCV[,1], tovar.P.cor, tulum$G.VCV[,1], tulum.G.cor, tulum$E.VCV[,1], tulum.E.cor, tulum$P.VCV[,1], tulum.P.cor)
	data.tab1.col2 <- list(tovar$G.VCV[,2], tovar$G.VCV[,4], tovar$E.VCV[,2], tovar$E.VCV[,4], tovar$P.VCV[,2], tovar$P.VCV[,4], tulum$G.VCV[,2], tulum$G.VCV[,4], tulum$E.VCV[,2], tulum$E.VCV[,4], tulum$P.VCV[,2], tulum$P.VCV[,4])
	
	# note that correlations should not be scaled
	tab1.col1 <- paste0(
		format(round(rep(c(mf,1), 6)*sapply(data.tab1.col1, mean), digits=dd),nsmall=dd), 
		" +/-",
		format(round(rep(c(mf,1), 6)*sapply(data.tab1.col1, sd), digits=dd),nsmall=dd),
		" (", 
		format(round(rep(c(mf,1), 6)*sapply(data.tab1.col1, function(x) HPDinterval(x)[1]), digits=dd), nsmall=dd),
		"; ", 
		format(round(rep(c(mf,1), 6)*sapply(data.tab1.col1, function(x) HPDinterval(x)[2]), digits=dd), nsmall=dd),
		")")
	
	tab1.col2 <- paste0(
		format(round(mf*sapply(data.tab1.col2, mean), digits=dd),nsmall=dd), 
		" +/-",
		format(round(mf*sapply(data.tab1.col2, sd), digits=dd),nsmall=dd),
		" (", 
		format(round(mf*sapply(data.tab1.col2, function(x) HPDinterval(x)[1]), digits=dd), nsmall=dd),
		"; ", 
		format(round(mf*sapply(data.tab1.col2, function(x) HPDinterval(x)[2]), digits=dd), nsmall=dd),
		")")
	
	table1 <- data.frame(' '=as.character(c("Tovar", "", "", "", "", "", "Tulum", "", "", "", "", "")), '  '=as.character(c("G-matrix", "", "E-matrix", "", "P-matrix", "", "G-matrix", "", "E-matrix", "", "P-matrix", "")), '   '=rep(c("GA","UBA"), 6), 'Gland area (GA)'=tab1.col1, 'Upper bract area (UBA)'=tab1.col2, check.names=FALSE)
	
	write.table(table1, file=if (!excluding.self) "table1a.txt" else "table1b.txt", quote=FALSE, sep="\t", row.names=FALSE)
}
