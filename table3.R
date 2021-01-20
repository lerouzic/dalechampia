#!/usr/bin/env Rscript

source("scripts/sraAutoreg-bivar.R")
source("scripts/summarypop.R")

makesradata <- function(dat, centering=c("raw","control", "updown")[1], exclude.control=TRUE) {
	# Formats the data from "summarypop" in a way that can be provided to the "sra" functions
	
	
	# Missing values in controls is a problem, replacing by the up/down average
	dat$mu.x[dat$Rep == "Control"] <- ifelse(is.na(dat$mu.x[dat$Rep == "Control"]), 0.5*dat$mu.x[dat$Rep == "Up"] + 0.5*dat$mu.x[dat$Rep == "Down"], dat$mu.x[dat$Rep == "Control"])
	dat$sig2.x[dat$Rep == "Control"] <- ifelse(is.na(dat$sig2.x[dat$Rep == "Control"]), 0.5*dat$sig2.x[dat$Rep == "Up"] + 0.5*dat$sig2.x[dat$Rep == "Down"], dat$sig2.x[dat$Rep == "Control"])
	dat$mu.y[dat$Rep == "Control"] <- ifelse(is.na(dat$mu.y[dat$Rep == "Control"]), 0.5*dat$mu.y[dat$Rep == "Up"] + 0.5*dat$mu.y[dat$Rep == "Down"], dat$mu.y[dat$Rep == "Control"])
	dat$sig2.y[dat$Rep == "Control"] <- ifelse(is.na(dat$sig2.y[dat$Rep == "Control"]), 0.5*dat$sig2.y[dat$Rep == "Up"] + 0.5*dat$sig2.y[dat$Rep == "Down"], dat$sig2.y[dat$Rep == "Control"])
	dat$sig.xy[dat$Rep == "Control"] <- ifelse(is.na(dat$sig.xy[dat$Rep == "Control"]), 0.5*dat$sig.xy[dat$Rep == "Up"] + 0.5*dat$sig.xy[dat$Rep == "Down"], dat$sig.xy[dat$Rep == "Control"])
	
	cor.x <- 0
	cor.y <- 0
	if (centering == "control") {
		cor.x <- dat$mu.x[dat$Rep=="Control"]
		cor.y <- dat$mu.y[dat$Rep=="Control"]
	}
	if (centering == "updown") {
		cor.x <- 0.5*dat$mu.x[dat$Rep=="Up"] + 0.5*dat$mu.x[dat$Rep=="Down"]
		cor.y <- 0.5*dat$mu.y[dat$Rep=="Up"] + 0.5*dat$mu.y[dat$Rep=="Down"]
	}
	
	if (exclude.control) 
		dat <- dat[dat$Rep != "Control",]
	
	if(is.null(dat$deltaV)) dat$deltaV <- rep(NA, nrow(dat))
	
	dat$S[is.na(dat$S)] <- 0
	
	ans <- data.frame(
		rep=factor(as.character(dat$Rep)),
		gen=dat$Gen, 
		phen.X.mean=dat$mu.x - cor.x,
		sel.X.mean=dat$mu.x+ifelse(is.na(dat$S), 0, dat$S) - cor.x,
		phen.X.var=dat$sig2.x,
		sel.X.var=dat$sig2s.x,
		phen.Y.mean=dat$mu.y - cor.y, 
		sel.Y.mean=dat$mu.y - cor.y,
		phen.Y.var=dat$sig2.y,
		sel.Y.var=dat$sig2.y, 
		phen.XY.covar=dat$sig.xy,
		N=dat$N)

	class(ans) <- c("sradata", class(ans))
	ans
}

maketable2 <- function(dat, CI.method=c("none", "vcov", "profile")[1], G0.boost=FALSE) {
	
	dataraw <- makesradata(dat, centering="raw")
	datacontrol <- makesradata(dat, centering="control")
	
	modsel <- function(data) {
		symnobul <- sraCstvar.bivar(data, Bulmer=FALSE, G0.boost=G0.boost)
		symbul   <- sraCstvar.bivar(data, Bulmer=TRUE, G0.boost=G0.boost)
		asymnobul<- sraCstvar.bivar.asym(data, Bulmer=FALSE, G0.boost=G0.boost)
		asymbul  <- sraCstvar.bivar.asym(data, Bulmer=TRUE, G0.boost=G0.boost)
		
		AICc   <- c(AICc.SRA(symnobul), AICc.SRA(symbul), AICc.SRA(asymnobul), NA, NA, AICc.SRA(asymbul), NA, NA)
		VA <- c(coef(symnobul)["logvarA0.A"], coef(symbul)["logvarA0.A"], NA, coef(asymnobul)["logvarA0.A.pos"], coef(asymnobul)["logvarA0.A.neg"], NA, coef(asymbul)["logvarA0.A.pos"], coef(asymbul)["logvarA0.A.neg"])
		covA  <- c(coef(symnobul)["covarA0.AB"], coef(symbul)["covarA0.AB"], NA, coef(asymnobul)["covarA0.AB.pos"], coef(asymnobul)["covarA0.AB.neg"], NA, coef(asymbul)["covarA0.AB.pos"], coef(asymbul)["covarA0.AB.neg"])
		
		if (CI.method=="vcov") {
			VA.se <- sqrt(c(vcov(symnobul)["logvarA0.A","logvarA0.A"], vcov(symbul)["logvarA0.A","logvarA0.A"], NA, vcov(asymnobul)["logvarA0.A.pos","logvarA0.A.pos"], vcov(asymnobul)["logvarA0.A.neg","logvarA0.A.neg"], NA, vcov(asymbul)["logvarA0.A.pos","logvarA0.A.pos"], vcov(asymbul)["logvarA0.A.neg","logvarA0.A.neg"]))
			covA.se <- sqrt(c(vcov(symnobul)["covarA0.AB","covarA0.AB"], vcov(symbul)["covarA0.AB","covarA0.AB"], NA, vcov(asymnobul)["covarA0.AB.pos","covarA0.AB.pos"], vcov(asymnobul)["covarA0.AB.neg","covarA0.AB.neg"], NA, vcov(asymbul)["covarA0.AB.pos","covarA0.AB.pos"], vcov(asymbul)["covarA0.AB.neg","covarA0.AB.neg"]))
			VA.min <- VA - 1.96*VA.se
			VA.max <- VA + 1.96*VA.se
			covA.min <- covA - 1.96*covA.se
			covA.max <- covA + 1.96*covA.se
		}
		data.frame(DeltaAICc=AICc-min(AICc, na.rm=TRUE), VA=VA, VA.min=VA.min, VA.max=VA.max, covA=covA, covA.min=covA.min, covA.max=covA.max)
	}

	transflogVA <- function(x) 100*exp(x)
	transfcovA <- function(x) 100*x
	
	modsel.raw <- modsel(dataraw)
	modsel.control <- modsel(datacontrol)
	
	tostring <- function(x, x.min=rep(NA,length(x)), x.max=rep(NA, length(x)), transf=identity, digits=2) {
		mapply(x, x.min, x.max, FUN=function(xx, xx.min, xx.max) {
			paste0(if(is.na(xx)) "" else as.character(round(transf(xx), digits=digits)), if(is.na(xx.min)) "" else paste0(" (", as.character(round(transf(xx.min), digits=digits)), "; ", if(transf(xx.max) > 100) "Inf" else as.character(round(transf(xx.max), digits=digits)), ")")) })
	}
	
	data.frame(
		Data=c("Raw", rep("", nrow(modsel.raw)-1), "Control-centered", rep("", nrow(modsel.control)-1)),
		Model=rep(c("Symmetric", "Symmetric + Bulmer", "Asymmetric", "", "", "Asymmetric + Bulmer", "", ""), 2),
		DeltaAICc=tostring(c(modsel.raw$DeltaAICc, modsel.control$DeltaAICc), digits=1), 
		'G.log.GAx100'=tostring(c(modsel.raw$VA, modsel.control$VA), x.min=c(modsel.raw$VA.min, modsel.control$VA.min), x.max=c(modsel.raw$VA.max, modsel.control$VA.max), transf=transflogVA, digits=2),
		'G.log.GA.log.UBAx100'=tostring(c(modsel.raw$covA, modsel.control$covA), x.min=c(modsel.raw$covA.min, modsel.control$covA.min), x.max=c(modsel.raw$covA.max, modsel.control$covA.max), transf=transfcovA, digits=2)
	)
}

cat("Tovar\n", file="table3.txt")
write.table(maketable2(summarypop("data/TovarData.txt"), CI.method="vcov"), file="table3.txt", append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
cat("Tulum\n", file="table3.txt", append=TRUE)
write.table(maketable2(summarypop("data/TulumData.txt"), CI.method="vcov", G0.boost=TRUE), file="table3.txt", append=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
