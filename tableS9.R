#!/usr/bin/env Rscript

source("scripts/plotpredict.R")
source("scripts/summarypop.R")


source("scripts/Gmatrices.R")


tovar.mat <- get.matrices("Tovar")
tulum.mat <- get.matrices("Tulum")

tovar <- summarypop("data/TovarData.txt")
tulum <- summarypop("data/TulumData.txt")

data.recenter.tovar.GA <- recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.x", normalization="raw")
data.recenter.tovar.UBA <- recenter(tovar, G=tovar.mat$G, Gv=tovar.mat$Gv, P=tovar.mat$P, N=64, Np=12, target="mu.y", normalization="raw")
data.recenter.tulum.GA <- recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.x", normalization="raw")
data.recenter.tulum.UBA <- recenter(tulum, G=tulum.mat$G, Gv=tulum.mat$Gv, P=tulum.mat$P, N=64, Np=12, target="mu.y", normalization="raw")

tab.tovar.GA <- data.frame('GA:V(sel)'=data.recenter.tovar.GA$va, 'GA:V(drift)'=data.recenter.tovar.GA$drift, 'GA:V(env)'=data.recenter.tovar.GA$env)
tab.tovar.UBA <- data.frame('UBA:V(sel)'=data.recenter.tovar.UBA$va, 'UBA:V(drift)'=data.recenter.tovar.UBA$drift, 'UBA:V(env)'=data.recenter.tovar.UBA$env)
tab.tulum.GA <- data.frame('GA:V(sel)'=data.recenter.tulum.GA$va, 'GA:V(drift)'=data.recenter.tulum.GA$drift, 'GA:V(env)'=data.recenter.tulum.GA$env)
tab.tulum.UBA <- data.frame('UBA:V(sel)'=data.recenter.tulum.UBA$va, 'UBA:V(drift)'=data.recenter.tulum.UBA$drift, 'UBA:V(env)'=data.recenter.tulum.UBA$env)

rownames(tab.tovar.GA) <- rownames(tab.tovar.UBA) <- rownames(tab.tulum.GA) <- rownames(tab.tulum.UBA) <- 0:4

ts9 <- "tableS9.txt"
mult <- 10000
digits <- 2

cat("Var x ", mult, "\nTovar\n", file=ts9)
write.table(format(round(mult*cbind(tab.tovar.GA, tab.tovar.UBA),digits=digits), nsmall=digits), sep="\t", quote=FALSE, file=ts9, append=TRUE)
cat("Tulum\n", file=ts9, append=TRUE)
write.table(format(round(mult*cbind(tab.tulum.GA, tab.tulum.UBA), digits=digits), nsmall=digits), sep="\t", quote=FALSE, file=ts9, append=TRUE)
