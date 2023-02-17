#!/usr/bin/env Rscript

source("../src/amodel.R")
source("../figures/common-colors.R")

# deleterious TEs - deleterious clusters

lty.max<- 3
lty.eq <- 2

N.sim <- 1000
rep.sim <- 20
Tmax <- 300

use.cache=TRUE
pdf("fig5.pdf", width=5.5, height=3.9, pointsize=7)

layout(rbind(1:2, 3:4))
par(mar=c(1,4,1,1), oma=c(3,0,0,0), cex=1)

model.default <- c(u=0.1, pi=0.03, s=0.01, k=1, sp=0.01, n0=1, p0=0)

model.par.s <- list(c(s=0, sp=0), c(s=0.01, sp=0.01), c(s=0.02, sp=0.02))
col.s       <- col[c("default", "s+", "s++")]

plot.model.dyn(model.default, model.par.s, what="n", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col.s)
axis(1, labels=FALSE)
axis(1, labels=FALSE)
text(10, 72, expression(hat(n)), col=col.s[1])
text(10, 10, expression(hat(n)), col=col.s[2])
text(10, 18, expression(hat(n)), col=col.s[3])


plot.model.dyn(model.default, model.par.s, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, xlab="", xaxt="n", col=col.s)
axis(1, labels=FALSE)
legend("bottomright", paste0("s=", sapply(model.par.s, "[", "s")), lty=1, col=col.s[seq_along(model.par.s)], bty="n")
text(10, 1.05, expression(hat(p)), col=col.s[1])
text(10, 0.71, expression(hat(p)), col=col.s[2])
text(10, 0.45, expression(hat(p)), col=col.s[3])


model.par.k <- list(c(k=1), c(k=2), c(k=5))
col.k       <- col[c("default","k2","k5")]

plot.model.dyn(model.default, model.par.k, what="n", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col.k)
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=1); par(xpd=FALSE)
text(10, 16, expression(hat(n)), col=col.k[1])
text(10, 30, expression(hat(n)), col=col.k[2])
text(10, 46, expression(hat(n)), col=col.k[3])

plot.model.dyn(model.default, model.par.k, what="p", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="topleft", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, col=col.k)
par(xpd=NA); mtext("Time (generations)", side=1, line=2.5, cex=1); par(xpd=FALSE)
text(10, 0.70, expression(hat(p)), col=col.k[1])
text(10, 0.46, expression(hat(p)), col=col.k[2])
text(10, 0.24, expression(hat(p)), col=col.k[3])

dev.off()

