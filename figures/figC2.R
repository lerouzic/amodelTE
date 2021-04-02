source("../src/amodel.R")

param.ref <- c(u=0.1, pi=0.03, k=1, s=0.01, sp=0)
init      <- c(n=1, p=0)

density <- 101

eqp <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p
}

s.expl <- seq(0.0001, 0.1, length=density)
u.expl <- seq(0.001, 1, length=density)
k.expl <- c(1, 2, 5)

s.eq.p <- lapply(k.expl, function(k) vapply(s.expl, function(s) {pp <- param.ref; pp["k"] <- k; pp["s"] <- s; eqp(pp, init) }, numeric(1)))
u.eq.p <- lapply(k.expl, function(k) vapply(u.expl, function(u) {pp <- param.ref; pp["k"] <- k; pp["u"] <- u; eqp(pp, init) }, numeric(1)))

pdf("figC2.pdf", width=8, height=4)
layout(t(1:2))
par(mar=c(5, 4.5, 1, 1))
plot(NULL, xlim=range(s.expl), ylim=c(0, 1), xlab="s", ylab=expression(hat(p)))
for (ki in seq_along(k.expl)) {
	lines(s.expl, s.eq.p[[ki]], lty=ki)
}
legend("topright", lty=seq_along(k.expl), col="darkgray", paste0("k = ", k.expl))


plot(NULL, xlim=range(u.expl), ylim=c(0, 1), xlab="u", ylab=expression(hat(p)))
for (ki in seq_along(k.expl)) {
	lines(u.expl, u.eq.p[[ki]], lty=ki)
}
dev.off()
