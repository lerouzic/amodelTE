source("../src/amodel.R")

library(viridis)

param.ref <- c(u=0.1, pi=0.03, k=2, s=0.01, sp=0.01)
init      <- c(n=1, p=0)
col       <- c(n="black", p="red", Re="black", Im="blue")
col[1:2] <- magma(3)[1:2]

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

ylim.n <- c(0,50)
ylim.p <- c(0,1)

pdf("fig5.pdf", width=9, height=6)
layout(matrix(1:6, ncol=3))
par(cex=1, mar=c(1,1,1,1), oma=c(5,4,0,0))

plot(NULL, xlim=range(s.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.n[[ki]], lty=ki, col=col["n"])
}
legend("topright", lty=seq_along(k.expl), col="darkgray", paste0("k = ", k.expl))
mtext(expression("Copy number ("*hat(n)*")"), 2, line=3, xpd=NA)
axis(1, labels=FALSE)
axis(2, xpd=NA)

plot(NULL, xlim=range(s.expl), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.p[[ki]], lty=ki, col=col["p"])
}
mtext(expression("Cluster frequency ("*hat(p)*")"), 2, line=3, xpd=NA)
mtext("Selection coefficient (s)", 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, xpd=NA)

plot(NULL, xlim=range(u.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.n[[ki]], lty=ki, col=col["n"])
}
axis(1, labels=FALSE)
axis(2, labels=FALSE)

plot(NULL, xlim=range(u.expl), ylim=ylim.p, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.p[[ki]], lty=ki, col=col["p"])
}
mtext("Transposition rate (u)", 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, labels=FALSE)

plot(NULL, xlim=range(pi.expl), ylim=ylim.n, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.n[[ki]], lty=ki, col=col["n"])
}
axis(1, labels=FALSE)
axis(2, labels=FALSE)

plot(NULL, xlim=range(pi.expl), ylim=ylim.p, xlab="", ylab="", xaxt="n", yaxt="n")
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.p[[ki]], lty=ki, col=col["p"])
}
mtext(expression("Cluster size ("*pi*")"), 1, line=3, xpd=NA)
axis(1, xpd=NA)
axis(2, labels=FALSE)

dev.off()
