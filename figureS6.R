#!/usr/bin/env Rscript

source("./scripts/plotpredict.R")

library(ellipse)

cols <- c(full="blue", exself="deeppink2")
transp <- 5

pdf("figureS6.pdf", width=10, height=5)

layout(t(1:2))

for (pop in c("Tovar", "Tulum")) {
	load(paste0("./data/", pop, "_MCMC.RData"))
	full   <- mod$VCV
	load(paste0("./data/", pop, "_MCMC_excludingself.RData"))
	exself <- mod$VCV
	
	genet.VCV <- grep(colnames(full), pattern="animal")
	full   <- full[,genet.VCV]
	exself <- exself[,genet.VCV]
	
	plot(NULL, xlim=c(-35,35), ylim=c(-35,35), xlab="log GA x 100", ylab="log UBA x 100", main=pop)
	
	for (i in 1:nrow(full))
		lines(ellipse(matrix(full[i,], ncol=2)), col=makeTransparent(cols["full"], transp))
	for (i in 1:nrow(exself))
		lines(ellipse(matrix(exself[i,], ncol=2)), col=makeTransparent(cols["exself"], transp))
	
	lines(ellipse(matrix((colMeans(full)), ncol=2)), col=cols["full"], lwd=3)
	lines(ellipse(matrix((colMeans(exself)), ncol=2)), col=cols["exself"], lwd=3)
	
	if (pop=="Tovar") 
		legend("topleft", lty=1, col=cols, legend=c("Full data", "W/o self"), lwd=3)
}

dev.off()
