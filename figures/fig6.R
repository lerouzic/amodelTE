#!/usr/bin/env Rscript

source("../src/amodel.R")
source("../figures/common-colors.R")


param.ref <- c(u=0.1, pi=0.03, k=1, s=0.01, sp=0.01)
init      <- c(n=1, p=0)

mtext.line <- 2.5

density <- 201

eqn <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$n
}

eqp <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p
}

s.expl <- seq(0.0, 0.1, length=density)
u.expl <- seq(0.001, 1, length=density)
pi.expl <-seq(0.0001,0.2, length=density) 
k.expl <- c(1, 2, 5)

s.eq.n <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; eqn(pp, init) }, numeric(1)))
s.eq.p <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; eqp(pp, init) }, numeric(1)))

u.eq.n <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqn(pp, init) }, numeric(1)))
u.eq.p <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqp(pp, init) }, numeric(1)))

pi.eq.n <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) {pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; eqn(pp, init) }, numeric(1)))
pi.eq.p <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) {pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; eqp(pp, init) }, numeric(1)))

ylim.n <- c(0,80)
ylim.p <- c(0,1)

col.k  <- col[c("default","k2","k5")]
col.s  <- col[c("default", "s+", "s++")]
col.u  <- col[c("u+","default","u-")]

pdf("fig6.pdf",  width=5.5, height=3.9, pointsize=7)
layout(matrix(1:6, ncol=3))
par(cex=1, mar=c(1,1,1,1), oma=c(5,4,0,0), cex=1)

plot(NULL, xlim=range(s.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.n[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(0, 0.01, 0.02), 
		y=c(	eqn({pp <- param.ref; pp["s"] <- pp["sp"] <- 0; pp}, init), 
				eqn({pp <- param.ref; pp["s"] <- pp["sp"]<- 0.01; pp}, init), 
				eqn({pp <- param.ref; pp["s"] <- pp["sp"]<- 0.02; pp}, init)), 
		pch=16, col=col.s)
legend("topright", lty=1, col=col.k, paste0("k = ", k.expl))
mtext(expression("Eq. copy number ("*hat(n)*")"), 2, line=3, xpd=NA)
axis(1, labels=FALSE)
axis(2, xpd=NA)

plot(NULL, xlim=range(s.expl), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.p[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(0, 0.01, 0.02), 
		y=c(	eqp({pp <- param.ref; pp["s"] <- pp["sp"] <- 0; pp}, init), 
				eqp({pp <- param.ref; pp["s"] <- pp["sp"]<- 0.01; pp}, init), 
				eqp({pp <- param.ref; pp["s"] <- pp["sp"]<- 0.02; pp}, init)), 
		pch=16, col=col.s)
mtext(expression("Eq. TE cluster insertion frequency ("*hat(p)*")"), 2, line=3, xpd=NA)
mtext("Selection coefficient (s)", 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, xpd=NA)

plot(NULL, xlim=range(u.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.n[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(param.ref["u"],param.ref["u"],param.ref["u"]), 
		y=c(	eqn({pp <- param.ref; pp["k"] <- 1; pp}, init), 
				eqn({pp <- param.ref; pp["k"] <- 2; pp}, init), 
				eqn({pp <- param.ref; pp["k"] <- 5; pp}, init)), 
		pch=16, col=col.k)
axis(1, labels=FALSE)
axis(2, labels=FALSE)

plot(NULL, xlim=range(u.expl), ylim=ylim.p, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.p[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(param.ref["u"],param.ref["u"],param.ref["u"]), 
		y=c(	eqp({pp <- param.ref; pp["k"] <- 1; pp}, init), 
				eqp({pp <- param.ref; pp["k"] <- 2; pp}, init), 
				eqp({pp <- param.ref; pp["k"] <- 5; pp}, init)), 
		pch=16, col=col.k)
mtext("Transposition rate (u)", 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, labels=FALSE)

plot(NULL, xlim=range(pi.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.n[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(param.ref["pi"],param.ref["pi"],param.ref["pi"]), 
		y=c(	eqn({pp <- param.ref; pp["k"] <- 1; pp}, init), 
				eqn({pp <- param.ref; pp["k"] <- 2; pp}, init), 
				eqn({pp <- param.ref; pp["k"] <- 5; pp}, init)), 
		pch=16, col=col.k)
axis(1, labels=FALSE)
axis(2, labels=FALSE)

plot(NULL, xlim=range(pi.expl), ylim=ylim.p, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.p[[ki]], lty=1, col=col.k[ki])
}
points(	x=c(param.ref["pi"],param.ref["pi"],param.ref["pi"]), 
		y=c(	eqp({pp <- param.ref; pp["k"] <- 1; pp}, init), 
				eqp({pp <- param.ref; pp["k"] <- 2; pp}, init), 
				eqp({pp <- param.ref; pp["k"] <- 5; pp}, init)), 
		pch=16, col=col.k)
mtext(expression("Total cluster size ("*pi*")"), 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, labels=FALSE)

dev.off()
