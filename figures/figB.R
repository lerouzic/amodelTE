#!/usr/bin/env Rscript

source("../src/amodel.R")

pdf("figB.pdf", width=8, height=12)

u.seq <- 10^seq(-3,0,length.out=21)
pi.seq <- seq(0.01,0.2,length.out=21)
n0.seq <- seq(0.01, 2, length.out=11)
k.seq <- 1:10

default.par <- list(u=0.1, pi=0.03, s=0, k=1, n0=1, selk=FALSE, p0=0)

get.eq <- function(default.par, focal.par, what="Eq") {
	pp <- default.par
	pp[names(focal.par)] <- focal.par
	do.call(pred.eq, pp)[[what]]$n
}

get.tmax <- function(default.par, focal.par, Tmax=1000, prox=0.9) {
	pp <- default.par
	pp[names(focal.par)] <- focal.par
	dd <- do.call(amodel, c(pp, list(Tmax=Tmax)))$n
	which(dd >= prox*max(dd))[1]
}

layout(rbind(1:2,3:4,5:6,7:8))

par(cex=1, mar=c(4,5,1,1))
plot(u.seq, sapply(u.seq, function(u) get.tmax(default.par, c(u=u)) ), xlab="u", ylab="Time to 90% max copies", type="l", log="xy")
plot(u.seq, sapply(u.seq, function(u) get.eq(default.par, c(u=u)) ), xlab="u", ylab="Copy number at equilibrium", type="l", log="x", ylim=c(0, 35))

plot(pi.seq, sapply(pi.seq, function(pi) get.tmax(default.par, c(pi=pi)) ), xlab=expression(pi), ylab="Time to 90% max copies", type="l")
plot(pi.seq, sapply(pi.seq, function(pi) get.eq(default.par, c(pi=pi)) ), xlab=expression(pi), ylab="Copy number at equilibrium", type="l")

plot(n0.seq, sapply(n0.seq, function(n0) get.tmax(default.par, c(n0=n0)) ), xlab=expression(n[0]), ylab="Time to 90% max copies", type="l")
plot(n0.seq, sapply(n0.seq, function(n0) get.eq(default.par, c(n0=n0)) ), xlab=expression(n[0]), ylab="Copy number at equilibrium", type="l", ylim=c(0, 35))

plot(k.seq, sapply(k.seq, function(k) get.tmax(default.par, c(k=k)) ), xlab="k", ylab="Time to 90% max copies", type="l", log="y")
plot(k.seq, sapply(k.seq, function(k) get.eq(default.par, c(k=k)) ), xlab="k", ylab="Copy number at equilibrium", type="l", log="", ylim=c(0, 400))

dev.off()
