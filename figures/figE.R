source("../src/amodel.R")

param.ref <- c(u=0.1, pi=0.03, k=2, s=0.01, sp=0.01)
init      <- c(n=1, p=0)
col       <- c(n="black", p="red", Re="black", Im="blue")

mtext.line <- 2.5

density <- 101

eqn <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$n
}

eqp <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p
}

s.expl <- seq(0.0001, 0.1, length=density)
u.expl <- seq(0.001, 1, length=density)
pi.expl <-seq(0.0001,0.2, length=density) 
k.expl <- c(1, 2, 5)

s.eq.n <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; eqn(pp, init) }, numeric(1)))
s.eq.p <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; eqp(pp, init) }, numeric(1)))

u.eq.n <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqn(pp, init) }, numeric(1)))
u.eq.p <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqp(pp, init) }, numeric(1)))

pi.eq.n <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) {pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; eqn(pp, init) }, numeric(1)))
pi.eq.p <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) {pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; eqp(pp, init) }, numeric(1)))

pdf("figE.pdf", width=15, height=5)
layout(rbind(1:3))
par(cex=1, mar=c(5,4,1,4))

plot(NULL, xlim=range(s.expl), ylim=c(0, 25), xlab="s", ylab="")
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.n[[ki]], lty=ki, col=col["n"])
}
par(new=TRUE)
plot(NULL, xlim=range(s.expl), ylim=c(0,1), xlab="", ylab="",  axes=FALSE)
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.p[[ki]], lty=ki, col=col["p"])
}
#~ axis(1, col=col["n"], col.axis=col["n"])
mtext(expression(hat(n)), side=2, line=mtext.line, col=col["n"], las=2)
axis(4, col=col["p"], col.axis=col["p"])
mtext(expression(hat(p)), side=4, line=mtext.line, col=col["p"], las=2)

legend("topright", lty=seq_along(k.expl), col="darkgray", paste0("k = ", k.expl))

plot(NULL, xlim=range(u.expl), ylim=c(0, 25), xlab="u", ylab="")
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.n[[ki]], lty=ki, col=col["n"])
}
par(new=TRUE)
plot(NULL, xlim=range(u.expl), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.p[[ki]], lty=ki, col=col["p"])
}
#~ axis(1, col=col["n"], col.axis=col["n"])
mtext(expression(hat(n)), side=2, line=mtext.line, col=col["n"], las=2)
axis(4, col=col["p"], col.axis=col["p"])
mtext(expression(hat(p)), side=4, line=mtext.line, col=col["p"], las=2)

plot(NULL, xlim=range(pi.expl), ylim=c(0, 25), xlab=expression(pi), ylab="")
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.n[[ki]], lty=ki, col=col["n"])
}
par(new=TRUE)
plot(NULL, xlim=range(pi.expl), ylim=c(0,1), xlab="", ylab="", axes=FALSE)
for (ki in seq_along(k.expl)) {
	lines(pi.expl, pi.eq.p[[ki]], lty=ki, col=col["p"])
}
#~ axis(1, col=col["n"], col.axis=col["n"])
mtext(expression(hat(n)), side=2, line=mtext.line, col=col["n"], las=2)
axis(4, col=col["p"], col.axis=col["p"])
mtext(expression(hat(p)), side=4, line=mtext.line, col=col["p"], las=2)

dev.off()
