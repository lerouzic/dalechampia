#!/usr/bin/env Rscript

# Generates population summary statistics (means/variances each generation)
# from the raw data. 

source("scripts/summarypop.R")

datadir <- "./data"

G0 <- 0

formatTable <- function(df, digits=3) {
	data.frame(
		'Rep.'=as.character(df$Rep),
		'Gen.'=as.character(df$Gen-1+G0),
		'mu.x'=format(round(df$mu.x, digits=digits), nsmall=digits),
		'sig2.x'=format(round(df$sig2.x, digits=digits), nsmall=digits),
		'S' = format(round(df$S, digits=digits), nsmall=digits),
		'deltaV'=format(round(df$deltaV, digits=digits), nsmall=digits),
		'mu.y'=format(round(df$mu.y, digits=digits), nsmall=digits),
		'sig2.y'=format(round(df$sig2.y, digits=digits), nsmall=digits),
		'sig.xy'=format(round(df$sig.xy, digits=digits), nsmall=digits),
		'N'=as.character(df$N),
		'SE.x'=format(round(df$se.x, digits=digits), nsmall=digits),
		'SE.y'=format(round(df$se.y, digits=digits), nsmall=digits)
	)
}

tovar.df <- summarypop(paste(datadir, "TovarData.txt", sep="/"))
tulum.df <- summarypop(paste(datadir, "TulumData.txt", sep="/"))

write.table(formatTable(tovar.df), file="tableS1.txt", row.names=FALSE, quote=FALSE, sep="\t", na="-")
write.table(formatTable(tulum.df), file="tableS2.txt", row.names=FALSE, quote=FALSE, sep="\t", na="-")
