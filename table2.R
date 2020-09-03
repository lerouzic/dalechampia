#!/usr/bin/env Rscript

source("scripts/plotpredict.R")
source("scripts/summarypop.R")

source("scripts/Gmatrices.R")

compute.all.vars <- function(data.recenter, mult=1) {
	
	ans <- list()
	ans$G       <- mult*data.recenter$va
	ans$D       <- mult*data.recenter$drift 
	ans$E       <-  mult*data.recenter$env
		
	ans$'Var#'  <- ans$G + ans$D/2 + ans$E/2
	ans$'%G#'   <- 100*ans$G / ans$'Var#'
	ans$'%D#'   <- 100*ans$D/2 / ans$'Var#'
	ans$'%E#'   <- 100*ans$E/2 / ans$'Var#'
	ans$'Resp.' <- (data.recenter$Up$pred - data.recenter$Down$pred)/2
	ans$'2SE#'  <- 2*sqrt(ans$'Var#'/mult)
	ans$'Rel.Err.#' <- sqrt(mult)*ans$'2SE#'/2/ans$'Resp.'
	
	ans$'Var$'  <- ans$G + 2*ans$D + 2*ans$E
	ans$'%G$'   <- 100*ans$G / ans$'Var$'
	ans$'%D$'   <- 100*ans$D*2 / ans$'Var$'
	ans$'%E$'   <- 100*ans$E*2 / ans$'Var$'
	ans$'2SE$'  <- 2*sqrt(ans$'Var$'/mult)
	ans$'Rel.Err.$' <- sqrt(mult)*ans$'2SE$'/2/ans$'Resp.'	
	
	as.data.frame(ans, check.names=FALSE)
}

format.mytable <- function(tt, digits.var=1, digits.resp=2) {
	for (cc in colnames(tt)) {
		if (grepl(cc, pattern="Resp") || grepl(cc, pattern="2SE"))
			tt[,cc] <- format(round(tt[,cc], digits=digits.resp), nsmall=digits.resp)
		else
			tt[,cc] <- format(round(tt[,cc], digits=digits.var), nsmall=digits.var)
	}
	tt
}


tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")

data.recenter.tovar.GA <- recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.x", normalization="raw")
data.recenter.tovar.UBA <- recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.y", normalization="raw")
data.recenter.tulum.GA <- recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.x", normalization="raw")
data.recenter.tulum.UBA <- recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.y", normalization="raw")

mult <- 10000

tab.tovar.GA  <- compute.all.vars(data.recenter.tovar.GA,  mult=mult)
tab.tovar.UBA <- compute.all.vars(data.recenter.tovar.UBA, mult=mult) 
tab.tulum.GA  <- compute.all.vars(data.recenter.tulum.GA,  mult=mult)
tab.tulum.UBA <- compute.all.vars(data.recenter.tulum.UBA, mult=mult)

rownames(tab.tovar.GA) <- rownames(tab.tovar.UBA) <- rownames(tab.tulum.GA) <- rownames(tab.tulum.UBA) <- 0:4


ts9 <- "table2.txt"
digits <- 1

cat("Var x ", mult, "\nTovar\nGA\n", file=ts9)
write.table(format.mytable(tab.tovar.GA), sep="\t", quote=FALSE, file=ts9, append=TRUE)
cat("UBA\n", file=ts9, append=TRUE)
write.table(format.mytable(tab.tovar.UBA), sep="\t", quote=FALSE, file=ts9, append=TRUE)
cat("Tulum\nGA\n", file=ts9, append=TRUE)
write.table(format.mytable(tab.tulum.GA), sep="\t", quote=FALSE, file=ts9, append=TRUE)
cat("UBA\n", file=ts9, append=TRUE)
write.table(format.mytable(tab.tulum.UBA), sep="\t", quote=FALSE, file=ts9, append=TRUE)
