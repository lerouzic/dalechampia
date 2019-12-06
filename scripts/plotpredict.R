
################ Helper functions for plotting things ##################

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

plot.ts.common <- function(pred, gen=seq_along(pred), verr=NULL, data, data.se=NULL, CI.factor=1.96, col.line="black", col.point=col.line, col.err=makeTransparent(col.line, alpha=50), pch=1, ...) {
	
		if (!is.null(verr)) {
			polygon(c(gen, rev(gen)), c(pred-CI.factor*sqrt(verr), rev(pred+CI.factor*sqrt(verr))), 
			border=NA, col=col.err)
		}
        lines(gen, pred, col=col.line, lwd=3)
        points(gen[!is.na(data)], data[!is.na(data)], pch=pch, lty=2, type="b", col=col.point)
        if(!is.null(data.se)) {
			arrows(x0=gen, y0=data-data.se, y1=data+data.se,
            length=0.1, angle=90, code=3, col=col.point)
        }
}

# Raw phenotypes, no normalization
plot.raw <- function(data, G, Gv, P, N, Np, target="mu.x", col.data=c(Control="darkgray", Up="black", Down="black"), pch=18, CI.factor=1.96,  ylab=target, xlab="Generations", ylim=NULL, verr.only=FALSE, ...) {
    stopifnot(target %in% c("mu.x", "mu.y"))
    
    if(is.null(ylim)) ylim <- 0.2*c(-1,1) + range(data[,target], na.rm=TRUE)
    
    beta <- data$S/data$sig2.x
    
    ggsel   <- if(target == "mu.x") G[1,1] else G[2,1]
    ggdrift <- if(target == "mu.x") G[1,1] else G[2,2]
    gv <- if(target == "mu.x") Gv[1,1] else Gv[2,1]
    ev <- if(target == "mu.x") P[1,1] else P[2,2]-G[2,2]
    
    gen <- seq(1, max(data$Gen))
    
    verr <- list(
		"Control"  = ev/64 + ggdrift * (1/64 + gen*1/64),
		"Up"       = ev/64 + ggdrift * (1/64 + gen*(1/16 - 1/128)) + gv * gen^2 *  mean(beta[data$Rep=="Up"], na.rm=TRUE)^2,
		"Down"     = ev/64 + ggdrift * (1/64 + gen*(1/16 - 1/128)) + gv * gen^2 *  mean(beta[data$Rep=="Down"], na.rm=TRUE)^2)
		
	if(verr.only) {
		return(data.frame(
			Line=rep(names(verr), each=length(gen)), 
			Generation=gen, 
			beta=c(beta[data$Rep=="Control"],beta[data$Rep=="Up"],beta[data$Rep=="Down"]),
			Ve=ev/64,
			Vdrift=c(ggdrift * (1/64 + gen*1/64), ggdrift * (1/64 + gen*(1/16 - 1/128)), ggdrift * (1/64 + gen*(1/16 - 1/128))),
			VvG=c( gen * 0, gv * gen^2 *  c(0,beta[data$Rep=="Up"][-5])^2, gv * gen^2 * c(0,beta[data$Rep=="Down"][-5])^2),
			Vtot=unlist(verr)
		))
	}
    

    plot(NULL, xlim=range(data$Gen), ylim=ylim, xlab=xlab, ylab=ylab, ...)
    # data
    for (i in seq_along(col.data)) {
        ii <- data$Rep==names(col.data)[i]
        
        bb <- cumsum(beta[ii])[-sum(ii)]
        if (all(is.na(bb))) bb <- rep(0, length(gen)) # Control lines: no beta
        else bb <- c(0, bb)

		myverr <- verr[[names(col.data)[i]]]
		pred <- data[min(which(ii)), target] + ggsel * bb 

		plot.ts.common(pred, gen=gen, verr=myverr, 
			data=data[ii,target], data.se=data[ii, if(target=="mu.x") "SEx" else "SEy"],
			CI.factor=CI.factor, col.line=col.data[i])
    }
}


# Control normalization
plot.control <- function(data, G, Gv, P, N, Np, target="mu.x", col.data=c(Control="darkgray", Up="black", Down="black"), pch=18, CI.factor=1.96,  ylab=target, xlab="Generations", ylim=NULL, verr.only=FALSE, ...) {
    stopifnot(target %in% c("mu.x", "mu.y"))
    
    if(is.null(ylim)) ylim <- 0.2*c(-1,1) + range(data[,target]-data[data$Rep=="Control",target], na.rm=TRUE)
    
    beta <- data$S/data$sig2.x
    
    ggsel   <- if(target == "mu.x") G[1,1] else G[2,1]
    ggdrift <- if(target == "mu.x") G[1,1] else G[2,2]
    gv <- if(target == "mu.x") Gv[1,1] else Gv[2,1]    
    ev <- if(target == "mu.x") P[1,1]-G[1,1] else P[2,2]-G[2,2]

    gen <- seq(1, max(data$Gen))
    
    verr <- list(
		"Control"  = rep(0, length(gen)),
		"Up"       = 2*ev/64 + ggdrift * (1/64 + gen*(1/16 - 1/128)) + gv * gen^2 *  mean(beta[data$Rep=="Up"], na.rm=TRUE)^2 + ggdrift * (1/64 + gen*1/64), 
		"Down"     = 2*ev/64 + ggdrift * (1/64 + gen*(1/16 - 1/128)) + gv * gen^2 *  mean(beta[data$Rep=="Down"], na.rm=TRUE)^2 + ggdrift * (1/64 + gen*1/64))
    
    if(verr.only) {
		return(data.frame(
			Line=rep(names(verr), each=length(gen)), 
			Generation=gen, 
			beta=c(beta[data$Rep=="Control"],beta[data$Rep=="Up"],beta[data$Rep=="Down"]),
			Ve=ev/64,
			Vdrift=c(ggdrift * (1/64 + gen*1/64), ggdrift * (1/64 + gen*(1/16 - 1/128)), ggdrift * (1/64 + gen*(1/16 - 1/128))),
			VvG=c( gen * 0, gv * gen^2 *  c(0,beta[data$Rep=="Up"][-5])^2, gv * gen^2 * c(0,beta[data$Rep=="Down"][-5])^2),
			Vtot=unlist(verr)
		))
	}

    plot(NULL, xlim=range(data$Gen), ylim=ylim, xlab=xlab, ylab=ylab, ...)
    # data
    for (i in seq_along(col.data)) {
        ii <- data$Rep==names(col.data)[i]
        
        bb <- cumsum(beta[ii])[-sum(ii)]
        if (all(is.na(bb))) bb <- rep(0, length(gen)) # Control lines: no beta
        else bb <- c(0, bb)

		myverr <- verr[[names(col.data)[i]]]
		pred <- ggsel * bb

		plot.ts.common(pred, gen=gen, verr=myverr, 
			data=data[ii,target]-data[data$Rep=="Control",target], data.se=data[ii, if(target=="mu.x") "SEx" else "SEy"],
			CI.factor=CI.factor, col.line=col.data[i])
    }
}


# Up-Down normalization
plot.updown <- function(data, G, Gv, P, N, Np, target="mu.x", col.data=c(Control="darkgray", Up="black", Down="black"), pch=18, CI.factor=1.96,  ylab=target, xlab="Generation", ylim=NULL, verr.only=FALSE, ...) {
    stopifnot(target %in% c("mu.x", "mu.y"))
    
    if(is.null(ylim)) ylim=0.2*c(-1,1) + range(data[,target]-0.5*data[data$Rep=="Up",target]-0.5*data[data$Rep=="Down",target], na.rm=TRUE)
    
    beta <- data$S/data$sig2.x
    
    ggsel   <- if(target == "mu.x") G[1,1] else G[2,1]
    ggdrift <- if(target == "mu.x") G[1,1] else G[2,2]
    gv <- if(target == "mu.x") Gv[1,1] else Gv[2,1]    
    ev <- if(target == "mu.x") P[1,1]-G[1,1] else P[2,2]-G[2,2]

    gen <- seq(1, max(data$Gen))
    
    verr.sel <- ev/64 + ggdrift * (1/64 + gen*(1/16 - 1/128))
    verr.c <- gv * gen^2 *  mean(abs(beta[data$Rep%in%c("Up","Down")]), na.rm=TRUE)^2
    verr.cont <- ev/64 + ggdrift * (1/64 + gen*1/64)
    verr <- list(
		"Control"  = verr.cont+verr.sel/2,
		"Up"       = verr.sel/2+verr.c, 
		"Down"     = verr.sel/2+verr.c)
    
	if(verr.only) {
		return(data.frame(
			Line=rep(names(verr), each=length(gen)), 
			Generation=gen, 
			beta=c(beta[data$Rep=="Control"],beta[data$Rep=="Up"],beta[data$Rep=="Down"]),
			Ve=ev/64,
			Vdrift=c(ggdrift * (1/64 + gen*1/64), ggdrift * (1/64 + gen*(1/16 - 1/128)), ggdrift * (1/64 + gen*(1/16 - 1/128))),
			VvG=c( gen * 0, gv * gen^2 *  c(0,beta[data$Rep=="Up"][-5])^2, gv * gen^2 * c(0,beta[data$Rep=="Down"][-5])^2),
			Vtot=unlist(verr)
		))
	}    
    

    plot(NULL, xlim=range(data$Gen), ylim=ylim, xlab=xlab, ylab=ylab, ...)
    # data
    for (i in seq_along(col.data)) {
        ii <- data$Rep==names(col.data)[i]
        
        bb <- c(0, cumsum(beta[ii])[-sum(ii)])
        bb[is.na(bb)] <- 0

		myverr <- verr[[names(col.data)[i]]]
		pred <- ggsel * bb

		plot.ts.common(pred, gen=gen, verr=myverr, 
			data=data[ii,target]-0.5*data[data$Rep=="Up",target]-0.5*data[data$Rep=="Down",target], data.se=data[ii, if(target=="mu.x") "SEx" else "SEy"],
			CI.factor=CI.factor, col.line=col.data[i])
    }
}
