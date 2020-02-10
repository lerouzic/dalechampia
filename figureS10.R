#!/usr/bin/env Rscript

source("scripts/sraAutoreg-bivar.R")
source("scripts/summarypop.R")
source("scripts/Gmatrices.R")

library(ellipse)


# Duplicated function with table 2
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

makefigS10 <- function(dat, what="var", G=NULL, col.G="black", col.raw="green", col.control="blue", xylim=c(-6.5,-4), xlab="log var(up)", ylab="log var(down)", legend=TRUE) {
	# Some duplcated code with maketable2? 
	
	dataraw <- makesradata(dat, centering="raw")
	datacontrol <- makesradata(dat, centering="control")
	
	if (what=="var") {
		param.sra <- "logvarA0.A"
		param.mcmc <- "logGA:logGA.animal"
	} else if (what == "cov") {
		param.sra <- "covar.A0.AB"
		param.mcmc <- "logGA:logUBA.animal"
	} else {stop()}
	
	addtofig <- function(model, param, ...) {
		name.up <- if(param %in% names(coef(model))) param else paste0(param, ".pos")
		name.down <- if(param %in% names(coef(model))) param else paste0(param, ".neg")
		
		param.up <- coef(model)[name.up]
		param.down <- coef(model)[name.down]
		var.up <- model$vcov[name.up,name.up]
		var.down <- model$vcov[name.down,name.down]
		cov.updown <- 0 # model$vcov[name.up,name.down]
		
		points(param.up, param.down, ...)
		lines(ellipse(cbind(c(var.up,cov.updown),c(cov.updown,var.down)), centre=c(param.up,param.down), level=0.95), ...)
	}

	plot(NULL, xlim=xylim, ylim=xylim, xlab=xlab, ylab=ylab, asp=1)

	addtofig(sraCstvar.bivar(dataraw, Bulmer=FALSE), param.sra, pch=1, col=col.raw)
	addtofig(sraCstvar.bivar(dataraw, Bulmer=TRUE), param.sra, pch=2, col=col.raw)
	addtofig(sraCstvar.bivar.asym(dataraw, Bulmer=FALSE), param.sra, pch=1, col=col.raw)
	addtofig(sraCstvar.bivar.asym(dataraw, Bulmer=TRUE), param.sra, pch=2, col=col.raw)
	
	addtofig(sraCstvar.bivar(datacontrol, Bulmer=FALSE), param.sra, pch=1, col=col.control)
	addtofig(sraCstvar.bivar(datacontrol, Bulmer=TRUE), param.sra, pch=2, col=col.control)
	addtofig(sraCstvar.bivar.asym(datacontrol, Bulmer=FALSE), param.sra, pch=1, col=col.control)
	addtofig(sraCstvar.bivar.asym(datacontrol, Bulmer=TRUE), param.sra, pch=2, col=col.control)
	
	if(!is.null(G)) {
		G.mean <- mean(log(G$G.VCV[,param.mcmc]))
		G.var <-  var(log(G$G.VCV[,param.mcmc]))
		points(G.mean, G.mean, pch=19, col=col.G)
		lines(ellipse(cbind(c(G.var,0),c(0,G.var)), centre=c(G.mean,G.mean), level=0.95), col=col.G)
	}
	abline(b=1, a=0, lty=3, col="gray")
	
	if (legend) {
		legend("bottomright", pch=c(19, 1, 2, 1, 2), col=c(col.G, col.raw, col.raw, col.control, col.control), legend=c("Diallel", "Raw, no Bulmer", "Raw, Bulmer", "Control, no Bulmer", "Control, Bulmer"))
	}
}

pdf("figureS10.pdf", width=5, height=10)
layout(1:2)
	makefigS10(summarypop("data/TovarData.txt"), G=get.matrices("Tovar"))
	title("Tovar")
	makefigS10(summarypop("data/TulumData.txt"), G=get.matrices("Tulum"))
	title("Tulum")
dev.off()
