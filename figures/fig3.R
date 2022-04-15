source("../src/amodel.R")

param.ref <- c(u=0.1, pi=0.03, k=1, s=0.01, sp=0)
init      <- c(n=1, p=0)

col.approx <- "gray"

density <- 101

eqp.a <- function(param, init) {
	c(p=pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p)
}

eqp.n <- function(param, init) {
	c(p=num.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p)
}

s.expl <- seq(0.0001, 0.1, length=density)
u.expl <- seq(0.001, 1, length=density)
k.expl <- c(1, 2, 5)

s.eq.p.a <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- s; eqp.a(pp, init) }, numeric(1)))
s.eq.p.n <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- s; eqp.n(pp, init) }, numeric(1)))

u.eq.p.a <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqp.a(pp, init) }, numeric(1)))
u.eq.p.n <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqp.n(pp, init) }, numeric(1)))




pdf("fig3.pdf", width=8, height=3.5)
layout(t(1:2))
par(mar=c(5, 4.5, 1, 1))
plot(NULL, xlim=range(s.expl), ylim=c(0, 1), xlab="s", ylab=expression(hat(p)))
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.p.a[[ki]], lty=ki, col=col.approx)
	lines(s.expl, s.eq.p.n[[ki]], lty=ki)	
}
legend("topright", lty=c(seq_along(k.expl), 1), col=c(rep("black", length(k.expl)), col.approx), c(paste0("k = ", k.expl), "Approx"), bty="n")


plot(NULL, xlim=range(u.expl), ylim=c(0, 1), xlab="u", ylab=expression(hat(p)))
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.p.a[[ki]], lty=ki, col=col.approx)
	lines(u.expl, u.eq.p.n[[ki]], lty=ki)
}
dev.off()
