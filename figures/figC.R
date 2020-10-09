#!/usr/bin/env Rscript

source("../src/amodel.R")


# neutral model
pdf("figC.pdf", width=8, height=4)

Tmax <- 150

lty.max<- 3
lty.eq <- 2

col <- 1:10

N.sim <- 500
rep.sim <- 6
Tmax <- 150

use.cache=TRUE

layout(rbind(1:2))
model.default <- c(u=0.1, pi=0.03, s=0.01, k=1, selk=FALSE, dom=TRUE, n0=1, p0=0)


###### Effect of transposition rate
model.par.s <- list(c(s=0), c(s=0.01), c(s=0.02))

plot.model.dyn(model.default, model.par.s, what="n", pred=TRUE, sim=TRUE, legend=TRUE, legend.pos="bottomright", Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)
plot.model.dyn(model.default, model.par.s, what="p", pred=TRUE, sim=TRUE, legend=FALSE, Tmax=Tmax, N=N.sim, rep=rep.sim, use.cache=use.cache)

dev.off()
