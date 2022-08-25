#!/usr/bin/env Rscript

source("../src/amodel.R")
source("../figures/common-colors.R")

library(vioplot)

model.neutral <- c(u=0.045, pi=0.03, s=0.00, k=2, sp=0.00, n0=1, p0=0)
model.dtnc    <- c(u=0.13, pi=0.03, s=0.01, k=2, sp=0.00, n0=1, p0=0)
model.dtdc    <- c(u=0.07, pi=0.03, s=0.01, k=2, sp=0.01, n0=1, p0=0)
model.regul   <- c(u=0.17, pi=0.00, s=0.01, k=2, sp=0.00, n0=1, p0=0, r=0.30)

use.cache <- TRUE

mycol=col[c('Neutral', 'Del TEs - neutral clusters', 'Del TEs - del clusters', 'Copy number regul')]
Nn <- c(10, 20, 50, 100, 200, 500, 1000)
rep <- 100

Tmax <- 100

ss.neutral <- lapply(Nn, function(N) {
	ans <- do.call(simmodel, c(as.list(model.neutral), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache, mean=FALSE)))
})

ss.dtnc <- lapply(Nn, function(N) {
	ans <- do.call(simmodel, c(as.list(model.dtnc), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache, mean=FALSE)))
})

ss.dtdc <- lapply(Nn, function(N) {
	ans <- do.call(simmodel, c(as.list(model.dtdc), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache, mean=FALSE)))
})

ss.regul <- lapply(Nn, function(N) {
	ans <- do.call(simmodel, c(as.list(model.regul), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache, mean=FALSE)))
})


names(ss.neutral) <- names(ss.dtnc) <- names(ss.dtdc) <- names(ss.regul) <- as.character(Nn)

ylim <- c(0, 80)

ss.all <- unlist(lapply(seq_along(ss.neutral), function(i) c(ss.neutral[i], ss.dtnc[i], ss.dtdc[i], ss.regul[i])), recursive=FALSE) 



pdf("fig6A.pdf", height=5, width=5)
	plot(Nn, sapply(ss.neutral, function(x) var(x$n[nrow(x$n),])/mean(x$n[nrow(x$n),])), ylim=c(0.02, 22), log="xy", xlab="Population size N", ylab=expression("Var("*bar(n)*")"/bar(n)), pch=19, col=mycol[1])
	points(Nn, sapply(ss.dtnc, function(x) var(x$n[nrow(x$n),])/mean(x$n[nrow(x$n),])), pch=19, col=mycol[2])
	points(Nn, sapply(ss.dtdc, function(x) var(x$n[nrow(x$n),])/mean(x$n[nrow(x$n),])), pch=19, col=mycol[3])
	points(Nn, sapply(ss.regul, function(x) var(x$n[nrow(x$n),])/mean(x$n[nrow(x$n),])), pch=19, col=mycol[4])
	
	vref <- var(ss.neutral[[1]]$n[nrow(ss.neutral[[1]]$n),])/mean(ss.neutral[[1]]$n[nrow(ss.neutral[[1]]$n),])
	
	lines(Nn, vref*Nn[1]/Nn, lty=2, col=mycol[1])
	legend("bottomleft", pch=19, col=mycol, legend=names(mycol), bty="n")
dev.off()


num.plot  <- 20
N.plot    <- 100
Tmax.plot <- 400

model.dtdc.plot  <- c(u=0.07, pi=0.03, s=0.01, k=2, sp=0.01, n0=1, p0=0)
model.regul.plot <- c(u=0.07, pi=0.00, s=0.01, k=2, sp=0.00, n0=1, p0=0, r=0.30)

ss.dtdc.plot  <-  do.call(simmodel, c(as.list(model.dtdc.plot),  list(N=N.plot, Tmax=Tmax.plot, rep=num.plot, use.cache=use.cache, mean=FALSE)))
ss.regul.plot <-  do.call(simmodel, c(as.list(model.regul.plot), list(N=N.plot, Tmax=Tmax.plot, rep=num.plot, use.cache=use.cache, mean=FALSE)))

pdf("fig6B.pdf", height=5, width=5)
	plot(NULL, xlim=c(0,Tmax.plot), ylim=c(0, 70), xlab="Generations", ylab=expression("Average copy number "*bar(n)))
	for (i in seq_len(num.plot)) {
		lines(ss.dtdc.plot$n[,i], lty=1, col=mycol[3])
		lines(ss.regul.plot$n[,i], lty=1, col=mycol[4])
	}
	legend("topleft", lty=1, col=mycol[3:4], legend=names(mycol)[3:4])
dev.off()


pdf("figSB2.pdf", height=4, width=12)
	# vioplot is not flexible, better to call plot before
	at <- 1:(4*length(ss.neutral)) + rep(0:(length(ss.neutral)-1), each=4)
	plot(NULL, xlim=range(at), ylim=ylim,  xlab="Population size (N)", ylab="Copy number (n)", xaxt="n")
	vioplot(lapply(ss.all, function(x) x$n[nrow(x$n),]), col=rep(mycol, length(ss.neutral)), at=at, add=TRUE)
	axis(1, at=2.5+(0:(length(ss.neutral)-1))*5, label=as.character(Nn))
	legend("topright", lty=1, lwd=6, col=mycol, legend=names(mycol))
	abline(h=mean(ss.neutral[[length(Nn)]]$n[nrow(ss.neutral[[length(Nn)]]$n),]), col="darkgray", lty=3)
dev.off()
