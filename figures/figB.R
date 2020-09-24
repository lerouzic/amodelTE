#!/usr/bin/env Rscript

source("amodel1.R")

pdf("figB.pdf", width=8, height=12)

u.seq <- 10^seq(-3,0,length.out=21)
pi.seq <- seq(0.01,0.2,length.out=21)
n0.seq <- seq(0.01, 2, length.out=11)
layout(rbind(1:2,3:4,5:6))
par(cex=1, mar=c(4,5,1,0))
plot(u.seq, sapply(u.seq, function(u) mod1.timetoeq(u, rho=0.03, thresh=0.90)$t ), xlab="u", ylab="Time to 90% regulation", type="l", log="xy")
plot(u.seq, sapply(u.seq, function(u) mod1.timetoeq(u, rho=0.03, thresh=0.90)$n ), xlab="u", ylab="Copy number at equilibrium", type="l", log="x", ylim=c(0, 35))

plot(pi.seq, sapply(pi.seq, function(pi) mod1.timetoeq(u=0.2, rho=pi, thresh=0.90)$t ), xlab=expression(pi), ylab="Time to 90% regulation", type="l", ylim=c(0,80))
plot(pi.seq, sapply(pi.seq, function(pi) mod1.timetoeq(u=0.2, rho=pi, thresh=0.90)$n ), xlab=expression(pi), ylab="Copy number at equilibrium", type="l")

plot(n0.seq, sapply(n0.seq, function(n0) mod1.timetoeq(u=0.2, rho=0.03, n0=n0, thresh=0.90)$t ), xlab=expression(n[0]), ylab="Time to 90% regulation", type="l", ylim=c(0,100))
plot(n0.seq, sapply(n0.seq, function(n0) mod1.timetoeq(u=0.2, rho=0.03, n0=n0, thresh=0.90)$n ), xlab=expression(n[0]), ylab="Copy number at equilibrium", type="l", ylim=c(0, 35))

dev.off()
