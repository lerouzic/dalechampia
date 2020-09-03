#!/usr/bin/env Rscript

# Reproduces table 3 from the manuscript
#
# The table reports regression slopes between log(UBA) as a response variable and log(GA) as a predictor variable, from a 
# mixed-effect model including the individual level as a random effect. 
#
# The table also reports correlation coefficients (and associated 95% confidence intervals) between log GA and log UBA

library(lme4)

dd.names <- c("P.Parent", paste0("F", 1:4, ".Up"), "P.Parent", paste0("F", 1:4, ".Down"), "P.Parent", paste0("F", 1:4, ".Ran"))

slope.df <- function(data) {
	
	slopes <- lapply(split(data, list(data$Gener, data$Line)), function(dd) {
			if (nrow(dd) == 0) {
				NA 
			} else { 
				coef(summary(lmer(logUBA ~ logGA + (1|ind), data=dd)))["logGA",c("Estimate", "Std. Error")] }
			})
	return(cbind(beta=sapply(slopes[dd.names], "[", 1), se=sapply(slopes[dd.names], "[", 2)))
}

cor.df <- function(data) {
	ddind <- do.call(rbind, by(data, data$ind, function(df) {
			data.frame(
				Line=df$Line[1], 
				Gener=df$Gener[1], 
				GA=mean(df$GA, na.rm=TRUE),
				logGA=mean(df$logGA, na.rm=TRUE), 
				UBA=mean(df$UBA, na.rm=TRUE),
				logUBA=mean(df$logUBA, na.rm=TRUE))}))
	cors <- lapply(split(ddind, list(ddind$Gener, ddind$Line)), function(dd) {
		if (nrow(dd) == 0) {
			NA
		} else {
			cc <- cor.test(log(dd$GA), log(dd$UBA))
			c(cc$estimate, cc$conf.int)
		}})
	return(cbind(r=sapply(cors[dd.names], "[", 1), CI.low=sapply(cors[dd.names], "[", 2), CI.high=sapply(cors[dd.names], "[", 3)))
}

format.table <- function(slopes, cors, digits=2) {
	
	bb <- paste0(format(round(slopes[,1], digits=digits), nsmall=digits), " (±", format(round(slopes[,2], digits=digits), nsmall=digits), ")")
	rr <- paste0(format(round(cors[,1], digits=digits), nsmall=digits), " (", format(round(cors[,2],digits=digits), nsmall=digits), "; ", format(round(cors[,3],digits=digits), nsmall=digits), ")")
	
	mm <- cbind(bb[1:5], rr[1:5], bb[6:10], rr[6:10], bb[11:15], rr[11:15])
	rownames(mm) <- as.character(1:5) # c("Parents", paste0("F", 1:4))
	mm
}

tovar.raw <- read.table("./data/TovarData.txt", header=TRUE)
tulum.raw <- read.table("./data/TulumData.txt", header=TRUE)

cat("\tUp-selected lines\t\tDown-selected lines\t\tControl\t\n", file="table3.txt")
cat("\tβ (±SE)\tr (95% CI)\tβ (±SE)\tr (95% CI)\tβ (±SE)\tr (95% CI)\n", file="table4.txt", append=TRUE)
cat("Tovar\n", file="table4.txt", append=TRUE)
write.table(format.table(slope.df(tovar.raw), cor.df(tovar.raw)), file="table4.txt", col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
cat("Tulum\n", file="table4.txt", append=TRUE)
write.table(format.table(slope.df(tulum.raw), cor.df(tulum.raw)), file="table4.txt", col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
