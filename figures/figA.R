#!/usr/bin/env Rscript

source("../src/amodel.R")


# neutral model
pdf("figA.pdf", width=8, height=8)

Tmax <- 150

lty.max<- 3
lty.eq <- 2

col <- 1:10

N.sim <- 10000
rep.sim <- 6
Tmax <- 150

use.cache=TRUE

layout(rbind(1:2, 3:4))
model.default <- c(u=0.1, pi=0.03, s=0, k=1, sp=0, dom=TRUE, n0=1, p0=0)


###### Effect of transposition rate
model.par.u <- list(c(u=0.2), c(u=0.1), c(u=0.05))

plot.model.dyn(model.default, model.par.u, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)
plot.model.dyn(model.default, model.par.u, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE)

model.par.pi <- list(c(pi=0.01), c(pi=0.03), c(pi=0.1))

plot.model.dyn(model.default, model.par.pi, what="n", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)
plot.model.dyn(model.default, model.par.pi, what="p", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE)

dev.off()
