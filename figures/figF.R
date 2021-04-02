source("../src/amodel.R")

param.ref <- c(u=0.1, pi=0.03, k=1, s=0.01, sp=0.01)
init      <- c(n=1, p=0)
col       <- c(n="black", p="red", Re="black", Im="blue")

density <- 101

jacob.fun <- jacob.dtdc


ev1 <- function(param, init) {
	jj <- jacob.fun(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])
	eigen(jj)$values[1]
}

eqp <- function(param, init) {
	pred.eq(u=param["u"], pi=param["pi"], k=param["k"], s=param["s"], sp=param["sp"], n0=init["n"], p0=init["p"])$Eq$p
}

s.expl <- seq(0.0001, 0.3, length=density)
u.expl <- seq(0.001, 1, length=density)
pi.expl <-seq(0.0001,0.2, length=density) 
k.expl <- c(1, 2, 5)

s.1D <- lapply(k.expl, function(k) vapply(s.expl, function(s) { pp <- param.ref; pp["k"] <- k; pp["s"] <- pp["sp"] <- s; ev1(pp, init) }, complex(1)))
u.1D <- lapply(k.expl, function(k) vapply(u.expl, function(u) { pp <- param.ref; pp["k"] <- k; pp["u"] <- u; ev1(pp, init) }, complex(1)))
pi.1D <- lapply(k.expl, function(k) vapply(pi.expl, function(pi) { pp <- param.ref; pp["k"] <- k; pp["pi"] <- pi; ev1(pp, init) }, complex(1)))

su.2D <- outer(s.expl, u.expl, FUN=function(ss, uu) mapply(ss, uu, FUN=function(s, u) { pp <- param.ref; pp[c("s","sp","u")] <- c(s, s, u); ev1(pp, init) }))
eqp.2D <- outer(s.expl, u.expl, FUN=function(ss, uu) mapply(ss, uu, FUN=function(s, u) { pp <- param.ref; pp[c("s","sp","u")] <- c(s, s, u); eqp(pp, init) }))

su.2D[eqp.2D <= 0] <- NA
su.2D[eqp.2D > 1] <- NA # never happens


pdf("figF.pdf", width=15, height=5)
layout(rbind(1:3))
par(cex=1)

plot(NULL, xlim=range(s.expl), ylim=c(-0.1, 0.05), xlab="s", ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(s.expl, Re(s.1D[[ki]]), lty=ki, col=col["Re"])
	lines(s.expl, Im(s.1D[[ki]]), lty=ki, col=col["Im"])
}
legend("bottomleft", lty=c(seq_along(k.expl), 1, 1), col=c(rep("darkgray", length(k.expl)) , col["Re"], col["Im"]), legend=c(paste0("k=",k.expl), "Real","Imaginary"))


plot(NULL, xlim=range(u.expl), ylim=c(-0.03, 0.03), xlab="u", ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(u.expl, Re(u.1D[[ki]]), lty=ki, col=col["Re"])
	lines(u.expl, Im(u.1D[[ki]]), lty=ki, col=col["Im"])
}


plot(NULL, xlim=range(pi.expl), ylim=c(-0.03, 0.03), xlab=expression(pi), ylab=expression("First Eigenvalue "*lambda[1]))
for (ki in seq_along(k.expl)) {
	lines(pi.expl, Re(pi.1D[[ki]]), lty=ki, col=col["Re"])
	lines(pi.expl, Im(pi.1D[[ki]]), lty=ki, col=col["Im"])
}

dev.off()


pdf("figF2.pdf", width=5, height=5)

image(x=u.expl, y=s.expl, z=Re(t(su.2D)), col = hcl.colors(1024, "viridis", rev = FALSE), xlab="u", ylab="s")
contour(x=u.expl, y=s.expl, z=Re(t(su.2D)), add=TRUE, labcex=0.8)
contour(x=u.expl, y=s.expl, z=Im(t(su.2D)), levels=0.0001, lwd=3, add=TRUE, col="red")
curve(x/(1+2*x), add=TRUE, lwd=3, col="darkviolet")

text(0.6, 0.13, "Stable focus", cex=1.6, col="white")
text(0.4, 0.205, "Stable eq", cex=1.6, col="white", srt=45)
text(0.2, 0.25, "No eq", cex=1.6)


dev.off()
