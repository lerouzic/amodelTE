#!/usr/bin/env Rscript

source("../src/amodel.R")

library("viridis")

# neutral model
lty.max<- 3
lty.eq <- 2

col <- magma(4)

#~ N.sim <- 10000
#~ rep.sim <- 1

N.sim <- 1000
rep.sim <- 20

Tmax <- 300

use.cache=TRUE
model.default <- c(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0)


pdf("fig1.pdf", width=8, height=6)
layout(rbind(1:2, 3:4))
par(mar=c(1,4,1,1), oma=c(3,0,0,0))

###### Effect of transposition rate
model.par.u <- list(c(u=0.2), c(u=0.1), c(u=0.05))

plot.model.dyn(model.default, model.par.u, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col[seq_along(model.par.u)])
axis(1, labels=FALSE)
text(10, 72, expression(hat(n)*"="), col=col[1]); text(28, 72, expression(hat(n)*"="), col=col[2]); text(42, 72, expression(hat(n)), col=col[3])

plot.model.dyn(model.default, model.par.u, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, xlab="", xaxt="n", col=col, ylim=c(0,1.1))
axis(1, labels=FALSE)
text(10, 1.05, expression(hat(p)*"="), col=col[1]); text(28, 1.05, expression(hat(p)*"="), col=col[2]); text(42, 1.05, expression(hat(p)), col=col[3])




model.par.k <- list(c(k=1), c(k=2), c(k=3))

plot.model.dyn(model.default, model.par.k, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col, ylim=c(0, 200))
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 58, expression(hat(n)), col=col[1]); text(10, 125, expression(hat(n)), col=col[2]); text(10, 190, expression(hat(n)), col=col[3])


plot.model.dyn(model.default, model.par.k, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, col=col, ylim=c(0,1.1))
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 1.05, expression(hat(p)*"="), col=col[1]); text(28, 1.05, expression(hat(p)*"="), col=col[2]); text(42, 1.05, expression(hat(p)), col=col[3])


dev.off()
