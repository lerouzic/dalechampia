#!/usr/bin/env Rscript

source("scripts/summarypop.R")

datadir <- "data/"

beta.from.df <- function(df) 
	df$sig.xy / df$sig2.x
	
r.from.df <- function(df)
	df$sig.xy / sqrt(df$sig2.x * df$sig2.y)

rCI.from.df <- function(df, probs=c(0.025, 0.975)) {
	# Based on Fisher z transformation
	r <- r.from.df(df)
	z <- 0.5*log((1+r)/(1-r))
	se <- 1/sqrt(df$N-3)
	sapply(qnorm(probs), function(pp) tanh(z+pp*se))
}

betase.from.df <- function(df) {
	r <- r.from.df(df)
	sqrt(1/(df$N-2)) * sqrt(df$sig2.x/df$sig2.y) * sqrt(1-r^2)
}

format.table <- function(df, digits=2) {
	beta <- beta.from.df(df)
	betase <- betase.from.df(df)
	r <- r.from.df(df)
	rCI <- rCI.from.df(df)
	
	bb <- paste0(format(round(beta, digits=digits), nsmall=digits), " (±", format(round(betase, digits=digits), nsmall=digits), ")")
	rr <- paste0(format(round(r, digits=digits), nsmall=digits), " (", format(round(rCI[,1],digits=digits), nsmall=digits), "; ", format(round(rCI[,2],digits=digits), nsmall=digits), ")")
	
	mm <- cbind(bb[1:5], rr[1:5], bb[6:10], rr[6:10], bb[11:15], rr[11:15])
	rownames(mm) <- c("Parents", paste0("F", 1:4))
	mm
}

tovar.df <- summarypop(paste(datadir, "TovarData.txt", sep="/"))
tulum.df <- summarypop(paste(datadir, "TulumData.txt", sep="/"))

cat("\tUp-selected lines\t\tDown-selected lines\t\tControl\t\n", file="table3.txt")
cat("\tβ (±SE)\tr (95% CI)\tβ (±SE)\tr (95% CI)\tβ (±SE)\tr (95% CI)\n", file="table3.txt", append=TRUE)
cat("Tovar\n", file="table3.txt", append=TRUE)
write.table(format.table(tovar.df), file="table3.txt", col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
cat("Tulum\n", file="table3.txt", append=TRUE)
write.table(format.table(tulum.df), file="table3.txt", col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
