#!/usr/bin/env Rscript

source("../src/amodel.R")


# neutral model

Tmax <- 150

lty.max<- 3
lty.eq <- 2

col <- 1:10

N.sim <- 10000
rep.sim <- 1

use.cache=TRUE

pdf("fig2.pdf", width=8, height=6)

layout(rbind(1:2, 3:4))
par(mar=c(1,4,1,1), oma=c(3,0,0,0))

model.default <- c(u=0.1, pi=0.03, s=0.01, k=1, sp=0, n0=1, p0=0)

model.par.s <- list(c(s=0), c(s=0.01), c(s=0.02))
plot.model.dyn(model.default, model.par.s, what="n", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n")
axis(1, labels=FALSE)
text(10, 36, expression(hat(n)))
text(10, 17, expression(n*"*"), col=2)


plot.model.dyn(model.default, model.par.s, what="p", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, xlab="", xaxt="n")
axis(1, labels=FALSE)
text(10, 1.05, expression(hat(p)))
text(10, 0.45, expression(p*"*"), col=3)


model.par.k <- list(c(k=1), c(k=2), c(k=5))
plot.model.dyn(model.default, model.par.k, what="n", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n")
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 14, expression(n*"*"))

plot.model.dyn(model.default, model.par.k, what="p", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE)
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 0.90, expression(hat(p)))
text(10, 0.70, expression(p*"*"))


dev.off()
