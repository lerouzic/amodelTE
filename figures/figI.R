#!/usr/bin/env Rscript

source("../src/amodel.R")

makeTransparent<-function(someColor, alpha=70)
{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# neutral model
lty.max<- 3
lty.eq <- 2

col <- 1:10

#~ N.sim <- 10000
#~ rep.sim <- 1

N.sim <- 1000
rep.sim <- 50

Tmax <- 10000

use.cache=TRUE
model.default <- c(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0)

uu <- 10^seq(-2,0,length.out=6)#11)
kk <- c(1, 2)

col.k <- setNames(1:length(kk), as.character(kk))

model.table <- expand.grid(u=uu, k=kk)
model.par <- lapply(as.data.frame(t(model.table)), function(sp) setNames(sp, nm=colnames(model.table)))

ylim <- c(0,100)
xlim <- range(uu)

pred.res <- lapply(model.par, function(mm) {
	pp <- model.default
	pp[names(mm)] <- mm
	do.call(pred.eq, as.list(pp))
})

sim.res <- lapply(model.par, function(mm) {
	pp <- model.default
	pp[names(mm)] <- mm
	do.call(simmodel, c(as.list(pp), list(mean=FALSE, N=N.sim, Tmax=Tmax, rep=rep.sim, use.cache=use.cache, simulator=simulator.default)))
})
final.copies <- sapply(sim.res, function(x) x$n[nrow(x$n),])

pdf("figI.pdf", width=5, height=4)
	plot(NULL, xlab="Transposition rate u", ylab=paste0("Copy number after ", Tmax, " generations"), xlim=xlim, ylim=ylim, log="x", col="gray")
	
	for (k in kk) {
		points(rep(model.table[model.table[,"k"] == k,"u"],each=rep.sim), final.copies[,model.table[,"k"] == k], col=makeTransparent(col.k[as.character(k)]),cex=0.5)
		points(model.table[model.table[,"k"] == k,"u"], colMeans(final.copies)[model.table[,"k"] == k], col=col.k[as.character(k)], pch=16)
		lines(model.table[model.table[,"k"] == k,"u"], sapply(pred.res[model.table[,"k"] == k], function(x) x$Eq$n), col=col.k[as.character(k)])
	}
dev.off()
