#!/usr/bin/env Rscript

source("../src/amodel.R")


# neutral model
pdf("figD.pdf", width=8, height=8)

Tmax <- 150

lty.max<- 3
lty.eq <- 2

col <- 1:10

N.sim <- 100
rep.sim <- 6
Tmax <- 500

use.cache=TRUE

layout(rbind(1:2,3:4))
model.default <- c(u=0.1, pi=0.03, s=0.01, k=1, selk=TRUE, n0=1, p0=0)

model.par.s <- list(c(s=0), c(s=0.01), c(s=0.02))
plot.model.dyn(model.default, model.par.s, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)
plot.model.dyn(model.default, model.par.s, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE)


model.par.k <- list(c(k=1), c(k=2), c(k=5))
plot.model.dyn(model.default, model.par.k, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)
plot.model.dyn(model.default, model.par.k, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=TRUE)

dev.off()
