# Computes the G matrices from the mcmc chains
# The script does not run the model on the diallel data, 
# this is not the topic from the current paper. 

datadir <- "./data"

# the way to collect the data is not entirely clean. The chain is 
# stored in two files for each population, that are stored in a "mod"
# variable in independent R sessions. 

get.VCV <- function(pop) {
	library(coda) # for as.mcmc

	load(paste0(datadir, "/", pop, "_logUBA_logGA_times100_1.RData"))
	VCV.part1 <- mod
	load(paste0(datadir, "/", pop, "_logUBA_logGA_times100_2.RData"))
	VCV.part2 <- mod
	return(as.mcmc(rbind(VCV.part1$VCV, VCV.part2$VCV)))
}


get.matrices <- function(pop, scale=c("data", "mcmc")[1]) {
	VCV <- get.VCV(pop)
	
	genet.VCV <- grep(colnames(VCV), pattern="animal")
	resid.VCV <- grep(colnames(VCV), pattern="units")
	indiv.VCV <- grep(colnames(VCV), pattern="ind")
	
	G.VCV <- VCV[,genet.VCV]
	P.VCV <- VCV[,genet.VCV]+VCV[,resid.VCV]+VCV[,indiv.VCV]
	
	# In the MCMC procedure, phenotypes are multiplied by 100, and variances by 100*100
	scale.factor <- 1
	if (scale == "data") scale.factor <- 1/10000
	
	ans <- list(
		G.VCV = G.VCV*scale.factor,
		P.VCV = P.VCV*scale.factor,
		G = matrix(apply(G.VCV*scale.factor, 2, median), ncol=2),
		P = matrix(apply(P.VCV*scale.factor, 2, median), ncol=2),
		Gv = matrix(apply(G.VCV*scale.factor, 2, var), ncol=2))
}
