
################ Plotting selection time series ##################

makeTransparent<-function(someColor, alpha=100)
{
	# From https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
	# Author: Nick Sabbe
	# Licence : CC-attribution-SA from the conditions of the website
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# Plots a single time series
plot.ts.common <- function(pred, gen=seq_along(pred), verr=NULL, data, data.se=NULL, CI.factor=1.96, col.line="black", col.point=col.line, col.err=makeTransparent(col.line, alpha=50), pch=1, ...) {
	# pred: the vector of predicted phenotype
	# gen : generation numbers (default: 1:5)
	# verr: vector of prediction error variance
	# data: data points
	# data.se: standard errors for the observed phenotypic means
	# CI.factor: factor by which sqrt(verr) should be multiplied to figure prediction intervals
	# ... graphical line, points, color options
	
	# Prediction error (shaded area)
	if (!is.null(verr)) {
		polygon(c(gen, rev(gen)), c(pred-CI.factor*sqrt(verr), rev(pred+CI.factor*sqrt(verr))), 
		border=NA, col=col.err)
	}
	
	# Prediction (plain line)
	lines(gen, pred, col=col.line, lwd=3)
	
	# Data points
	points(gen[!is.na(data)], data[!is.na(data)], pch=pch, lty=2, type="b", col=col.point)
	
	# Data error (error bars)
	if(!is.null(data.se) && data.se[1] > 0) {
		arrows(x0=gen, y0=data-data.se, y1=data+data.se,
	length=0.1, angle=90, code=3, col=col.point)
	}
}

# Estimate means and variances in raw and centered datasets
# Three ways to analyse the data:
# * raw:               phenotype as observed experimentally
# * control-centered:  phenotype centered on the control line
# * updown-centered:   phenotype centered on the average between up and down selection lines (symmetric response)
recenter <- function(data, G, Gv, P, N, Np, target="mu.x", normalization=c("raw", "control", "updown")[1]) {
	# data:          summary statistics (from summarypop script)
	# G:             G matrix
	# Gv:            Estimated variance of each element of the G matrix
	# P:             P matrix
	# N:             Number of individuals measured
	# Np:            Number of selected parents
	# target:        Either "mu.x" (selected trait) or "mu.y" (correlated trait)
	
	stopifnot(target %in% c("mu.x", "mu.y"))
	stopifnot(normalization %in% c("raw","control","updown"))
	
	se.target <- if(target=="mu.x") "se.x" else "se.y"

	# The different variances to account for depend on the target
	Gsel   <- if(target == "mu.x") G[1,1] else G[2,1]
	Gdrift <- if(target == "mu.x") G[1,1] else G[2,2]
	Gerr <- if(target == "mu.x") Gv[1,1] else Gv[2,1]
	Eerr <- if(target == "mu.x") P[1,1]-G[1,1] else P[2,2]-G[2,2]
	
	# Selection gradient (always on the selected trait, the gradient on the correlated trait is 0 by definition)
	beta <- data$S/data$sig2.x
		# mean gradients to compute prediction errors: a bit sloppy, but this should not really matter
	mean.beta   <- mean(abs(beta[data$Rep=="Up" | data$Rep == "Down"]), na.rm=TRUE)
	
	gen <- seq(1, max(data$Gen))
	

	pred.control <- 
		if (normalization == "raw") { rep(data[data$Rep=="Control" & data$Gen == 1, target], length(gen)) }
		else if (normalization == "control") { rep(0, length(gen) ) }
		else if (normalization == "updown") { c(0, cumsum(0.5*(beta[data$Rep=="Down"] + beta[data$Rep=="Up"]))[-length(gen)]*Gsel) }
	pred.up <- 
		if (normalization == "raw") { data[data$Rep=="Up" & data$Gen == 1, target] + c(0, cumsum(beta[data$Rep=="Up"])[-length(gen)]*Gsel) }
		else if (normalization == "control") { c(0, cumsum(beta[data$Rep=="Up"])[-length(gen)]*Gsel) }
		else if (normalization == "updown") { c(0, cumsum(0.5*(beta[data$Rep=="Up"] - beta[data$Rep=="Down"]))[-length(gen)]*Gsel) }
	pred.down <- 
		if (normalization == "raw") { data[data$Rep=="Down" & data$Gen == 1, target] + c(0, cumsum(beta[data$Rep=="Down"])[-length(gen)]*Gsel) }
		else if (normalization == "control") { c(0, cumsum(beta[data$Rep=="Down"])[-length(gen)]*Gsel) }
		else if (normalization == "updown") { c(0, cumsum(0.5*(beta[data$Rep=="Down"] - beta[data$Rep=="Up"]))[-length(gen)]*Gsel) }

	# Calculation of prediction variances
	# Three terms : environmental sampling error + drift + error in the additive variance estimate
	# Environmental sampling error : at the current generation, always Ve/N
	# Drift : cumulative over generations, + Va / N every generation (+ Va/N for generation 0)
	#         in selected lines, there are only Np parents, but the offspring number is normalized
	#         theory shows that the increase in variance is 1/Np - 1/2N every generation
	# Error on the estimate of Va : 
	#         this error cancels for drift (assuming overestimation is as likely as underestimation)
	#         but not for selected lines. A term in Var(Va) * t^2 * beta^2 needs to be considered.
	#         this error cumulates (quadratically) ver generations
	
	verr.drift.control <- Eerr/N + Gdrift*(1/N + (gen-1)*(1/N))
	verr.drift.sel     <- Eerr/N + Gdrift*(1/N + (gen-1)*(1/Np - 1/(2*N)))
	verr.va.sel        <- Gerr * (gen-1)^2 * mean.beta^2
	
	# If the data is control or up-down centered, the error variance is redistributed:
	verr.control <-
		if (normalization == "raw") { verr.drift.control }
		else if (normalization == "control") { rep(0, length(gen)) }
		else if (normalization == "updown") { verr.drift.control + verr.drift.sel/2 }
	verr.sel     <-
		if (normalization == "raw") { verr.drift.sel + verr.va.sel }
		else if (normalization == "control") { verr.drift.sel + verr.va.sel + verr.drift.control }
		else if (normalization == "updown") { verr.drift.sel/2 + verr.va.sel }

	phen.control <-
		if (normalization == "raw") { data[data$Rep=="Control",target] }
		else if (normalization == "control") { rep(0, length(gen)) }
		else if (normalization == "updown") { data[data$Rep=="Control",target] - 0.5*data[data$Rep=="Up",target] - 0.5*data[data$Rep=="Down",target] }
	phen.up <- 
		if (normalization == "raw") { data[data$Rep=="Up",target] }
		else if (normalization == "control") { data[data$Rep=="Up",target] - data[data$Rep=="Control",target] }
		else if (normalization == "updown") { 0.5*data[data$Rep=="Up",target] - 0.5*data[data$Rep=="Down",target] }
	phen.down <- 
		if (normalization == "raw") { data[data$Rep=="Down",target] }
		else if (normalization == "control") { data[data$Rep=="Down",target] - data[data$Rep=="Control",target] }
		else if (normalization == "updown") { 0.5*data[data$Rep=="Down",target] - 0.5*data[data$Rep=="Up",target] }

	se.control <- 
		if (normalization == "raw") { data[data$Rep=="Control",se.target] }
		else if (normalization == "control") { rep(0, length(gen)) }
		else if (normalization == "updown") { sqrt(data[data$Rep=="Control",se.target]^2 + (1/4)*(data[data$Rep=="Up",se.target]^2 + data[data$Rep=="Down",se.target]^2)) }
	se.up <- 
		if (normalization == "raw") { data[data$Rep=="Up",se.target] }
		else if (normalization == "control") { sqrt(data[data$Rep=="Up",se.target]^2 + data[data$Rep=="Control",se.target]^2) }
		else if (normalization == "updown") { sqrt((1/4)*(data[data$Rep=="Up",se.target]^2 + data[data$Rep=="Down",se.target]^2)) }
	se.down <- 
		if (normalization == "raw") { data[data$Rep=="Down",se.target] }
		else if (normalization == "control") { sqrt(data[data$Rep=="Down",se.target]^2 + data[data$Rep=="Control",se.target]^2) }
		else if (normalization == "updown") { sqrt((1/4)*(data[data$Rep=="Up",se.target]^2 + data[data$Rep=="Down",se.target]^2)) }

	
	return(list(
		Control = data.frame(gen=gen, pred=pred.control, phen=phen.control, se=se.control, verr=verr.control),
		Up      = data.frame(gen=gen, pred=pred.up,      phen=phen.up,      se=se.up,      verr=verr.sel),
		Down    = data.frame(gen=gen, pred=pred.down,    phen=phen.down,    se=se.down,    verr=verr.sel)))
}


# Call the plot routine on the recentered data. 
plot.data.recenter <- function(data.recenter, col.data=c(Control="gray50", Up="black", Down="black"), pch=18, CI.factor=1.96,  ylab="Phenotype", xlab="Generations", ylim=NULL, ...) {

	if(is.null(ylim)) 
		ylim <- 0.2*c(-1,1) + range(do.call(c, lapply(data.recenter, function(x) x$phen)), na.rm=TRUE)

	plot(NULL, xlim=range(data.recenter$Control$gen), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	for (ll in names(data.recenter)) {
		plot.ts.common(data.recenter[[ll]]$pred, data.recenter[[ll]]$gen, data.recenter[[ll]]$verr, data.recenter[[ll]]$phen, data.recenter[[ll]]$se, CI.factor=CI.factor, col.line=col.data[ll])
	}
}
