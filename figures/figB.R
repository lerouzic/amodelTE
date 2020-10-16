#!/usr/bin/env Rscript

source("../src/amodel.R")


log.par <- list(u=TRUE, pi=FALSE, k=FALSE, s=TRUE)
lo <- 41 # number of points

seqs <- list(
	u  = if (log.par$u) 10^seq(log10(0.01), log10(0.5), length.out=lo) else seq(0.01,0.5,length.out=lo),
	pi = seq(0.01,0.2,length.out=lo),
	k  = 1:5,
	s  = if (log.par$s) 10^seq(log10(0.001), log10(0.1), length.out=lo) else seq(0, 0.1, length.out=lo))

default.par <- list(
	neutral = list(u=0.1, pi=0.03, s=0, k=1, n0=1, sp=0, p0=0),
	selTE   = list(u=0.1, pi=0.03, s=0.01, k=1, n0=1, sp=0, p0=0),
	selcl   = list(u=0.1, pi=0.03, s=0.01, k=1, n0=1, sp=0.01, p0=0))

col <- c(neutral="black", selTE="blue", selcl="red")
caption <- c(u="Transposition rate (u)", pi=expression("Cluster size ("*pi*")"), k="Number of clusters (k)", s="Selection coeffcient (s)")
 

get.eq <- function(default.par, focal.par, what="n", what.eq="Eq") {
	pp <- default.par
	pp[names(focal.par)] <- focal.par
	if (names(focal.par) == "s" && default.par[["sp"]] != 0)
		pp["sp"] <- focal.par
	do.call(pred.eq, pp)[[what.eq]][[what]]
}

get.tmax <- function(default.par, focal.par, Tmax=10000, prox=0.9) {
	pp <- default.par
	pp[names(focal.par)] <- focal.par
	if (names(focal.par) == "s" && default.par[["sp"]] != 0)
		pp["sp"] <- focal.par
	dd <- do.call(amodel, c(pp, list(Tmax=Tmax)))$n
	which(dd >= prox*max(dd))[1]
}

get.nmax <- function(default.par, focal.par, Tmax=10000, prox=0.9) {
	pp <- default.par
	pp[names(focal.par)] <- focal.par
	if (names(focal.par) == "s" && default.par[["sp"]] != 0)
		pp["sp"] <- focal.par
	dd <- do.call(amodel, c(pp, list(Tmax=Tmax)))$n
	max(dd)
}


pdf("figB.pdf", width=12, height=12)
layout(matrix(1:(4*length(seqs)), byrow=FALSE, ncol=4))

par(cex=1, mar=c(0.5,0.5,0.5,0.5), oma=c(5,4,0,0), lwd=2)

for (ppar in names(seqs)) {
	plot(NULL, xlim=range(seqs[[ppar]]), ylim=c(0,500), xlab="", ylab="", log=if(log.par[[ppar]]) "x" else "", xaxt="n", yaxt="n")
	if (ppar == names(seqs)[1]) 
		mtext("Time to 90% max copies", 2, line=2.5, xpd=NA)
	axis(1, labels=FALSE)
	axis(2, labels=(ppar == names(seqs)[1]))
	for (pp in names(col)) 
		lines(seqs[[ppar]], sapply(seqs[[ppar]], function(x) get.tmax(default.par[[pp]], setNames(x, ppar))), col=col[[pp]])

	if (ppar == names(seqs)[2])
		legend("topright", lty=1, col=col, legend=c("Neutral", "Del TEs", "Del TEs & Clusters"), bty="n")
	
	plot(NULL, xlim=range(seqs[[ppar]]), ylim=c(0,1), xlab="", ylab="", log=if(log.par[[ppar]]) "x" else "", xaxt="n", yaxt="n")
	if (ppar == names(seqs)[1]) 
		mtext(expression("Equilibrium custer frequency ("*hat(p)*")"), 2, line=2.5, xpd=NA)
	axis(1, labels=FALSE)
	axis(2, labels=(ppar == names(seqs)[1]))
	for (pp in names(col)) 
		lines(seqs[[ppar]], sapply(seqs[[ppar]], function(x) get.eq(default.par[[pp]], setNames(x, ppar), what="p")), col=col[[pp]])
		
	plot(NULL, xlim=range(seqs[[ppar]]), ylim=c(0,100), xlab="", ylab="", log=if(log.par[[ppar]]) "x" else "", xaxt="n", yaxt="n")
	if (ppar == names(seqs)[1]) 
		mtext(expression("Max copy number ("*n^"*"*")"), 2, line=2.5, xpd=NA)
	axis(1, labels=FALSE)
	axis(2, labels=(ppar == names(seqs)[1]))	
	for (pp in names(col)) 
		lines(seqs[[ppar]], sapply(seqs[[ppar]], function(x) get.nmax(default.par[[pp]], setNames(x, ppar))), col=col[[pp]])
	
	plot(NULL, xlim=range(seqs[[ppar]]), ylim=c(0,100), xlab="", ylab="", log=if(log.par[[ppar]]) "x" else "", xaxt="n", yaxt="n")
	mtext(caption[ppar], 1, line=2.5, xpd=NA)
	axis(1)
	axis(2, labels=(ppar == names(seqs)[1]))
	if (ppar == names(seqs)[1]) 
		mtext(expression("Equilibrium copy number ("*hat(n)*")"), 2, line=2.5, xpd=NA)
	for (pp in names(col)) 
		lines(seqs[[ppar]], sapply(seqs[[ppar]], function(x) get.eq(default.par[[pp]], setNames(x, ppar))), col=col[[pp]])	
}


dev.off()
