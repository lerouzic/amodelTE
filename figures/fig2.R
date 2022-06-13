#!/usr/bin/env Rscript

source("../src/amodel.R")

library("viridis")

# neutral model

lty.max<- 3
lty.eq <- 2

col <- magma(4)

N.sim <- 1000
rep.sim <- 20

use.cache=TRUE

Tmax <- 150

pdf("fig2.pdf", width=8, height=6)

layout(rbind(1:2, 3:4))
par(mar=c(1,4,1,1), oma=c(3,0,0,0))

model.default <- c(u=0.1, pi=0.03, s=0.01, k=1, sp=0, n0=1, p0=0)

model.par.s <- list(c(s=0), c(s=0.01), c(s=0.02))
plot.model.dyn(model.default, model.par.s, what="n", pred=TRUE, sim=TRUE, max=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col)
axis(1, labels=FALSE)
text(10, 72, expression(hat(n)), col=col[1])
text(10, 32, expression(n*"*"), col=col[2])
text(10, 21, expression(n*"*"), col=col[3])
text(70, 3, expression(hat(n)*"="), col=col[2])
text(77, 3, expression(hat(n)), col=col[3])


plot.model.dyn(model.default, model.par.s, what="p", pred=TRUE, sim=TRUE, max=FALSE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, xlab="", xaxt="n", col=col, ylim=c(0,1.1))
axis(1, labels=FALSE)
text(10, 1.05, expression(hat(p)), col=col[1])
text(10, 0.88, expression(hat(p)), col=col[2])
text(10, 0.73, expression(hat(p)), col=col[3])



model.par.k <- list(c(k=1), c(k=2), c(k=5))
plot.model.dyn(model.default, model.par.k, what="n", pred=TRUE, sim=TRUE, max=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache, xlab="", xaxt="n", col=col)
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 28, expression(n*"*"), col=col[1])
text(10, 33.5, expression(n*"*"), col=col[2])
text(10, 41, expression(n*"*"), col=col[3])
text(70, 3, expression(hat(n)*"="), col=col[1])
text(80, 3, expression(hat(n)*"="), col=col[2])
text(88, 3, expression(hat(n)), col=col[3])

plot.model.dyn(model.default, model.par.k, what="p", pred=TRUE, sim=TRUE, max=FALSE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE, col=col)
par(xpd=NA); axis(1); mtext("Time (generations)", side=1, line=2.5, cex=0.8); par(xpd=FALSE)
text(10, 0.90, expression(hat(p)), col=col[1])
text(10, 0.60, expression(hat(p)), col=col[2])
text(10, 0.31, expression(hat(p)), col=col[3])

dev.off()
