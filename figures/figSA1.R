#!/usr/bin/env Rscript

source("../src/amodel.R")
source("../figures/common-colors.R")

param.ref <- c(u=0.1, pi=0.03, k=1, s=0.01, sp=0.01)
init      <- c(n=1, p=0)
col.k     <- col[c("default","k2","k5")]

density <- 101

jacob.fun <- jacob.dtdc


ev1 <- function(param, init) {
	jj <- jacob.fun(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])
	eigen(jj)$values[1]
}

eqp <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p
}

s.expl <- seq(0.0001, 0.08, length=density)
u.expl <- seq(0.001, 1, length=density)
pi.expl <-seq(0.0001,0.2, length=density) 
k.expl <- c(1, 2, 5)

s.1D <- lapply(k.expl, function(k) vapply(s.expl, function(s) { pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; ev1(pp, init) }, complex(1)))
u.1D <- lapply(k.expl, function(k) vapply(u.expl, function(u) { pp <- param.ref; pp["k"] <- k; pp["u"] <- u; ev1(pp, init) }, complex(1)))
pi.1D <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) { pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; ev1(pp, init) }, complex(1)))



pdf("figSA1.pdf", height=3.5, width=7.5, pointsize=7)
layout(rbind(1:3))
par(cex=1)

plot(NULL, xlim=range(s.expl), ylim=c(-0.1, 0.05), xlab="s", ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(s.expl, Re(s.1D[[ki]]), lty=1, col=col.k[ki])
	lines(s.expl, Im(s.1D[[ki]]), lty=2, col=col.k[ki])
}
legend("bottomleft", lty=c(1, 1, 1, 1, 2), col=c(col.k, rep("darkgray", 2)), legend=c(paste0("k=",k.expl), "Real","Imaginary"))


plot(NULL, xlim=range(u.expl), ylim=c(-0.03, 0.03), xlab="u", ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(u.expl, Re(u.1D[[ki]]), lty=1, col=col.k[ki])
	lines(u.expl, Im(u.1D[[ki]]), lty=2, col=col.k[ki])
}


plot(NULL, xlim=range(pi.expl), ylim=c(-0.03, 0.03), xlab=expression(pi), ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(pi.expl, Re(pi.1D[[ki]]), lty=1, col=col.k[ki])
	lines(pi.expl, Im(pi.1D[[ki]]), lty=2, col=col.k[ki])
}

dev.off()
