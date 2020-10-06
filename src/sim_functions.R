
library(parallel)

global.default <- list(
	mc.cores          = detectCores() - 1,
	mc.cores.internal = NA,
	nb.loci           = 1000,
	neutral.loci      = numeric(0), 
	piRNA.loci        = 1:30,
	fitness.FUN       = function(n.sel)    exp(-n.sel*0.01),
	regulation.FUN    = function(n.piRNA)  if (n.piRNA == 0) 1 else 0,
	u                 = 0.1,
	G                 = 100,
	N                 = 1000,
	rep               = 10,
	init.TE.ind       = 1
)

mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.cores=1) 
{
    FUN <- match.fun(FUN)
    answer <- mclapply(X = X, FUN = FUN, ..., mc.cores=mc.cores)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!isFALSE(simplify) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

n.sel.ind <- function(ind, global) {
	if (is.null(ind$n.sel))
		sum(! c(which(ind$gam1), which(ind$gam2)) %in% global$neutral.loci)
	else
		ind$n.sel
}

n.sel.pop <- function(pop, global) {
	mcsapply(pop, n.sel.ind, global=global, mc.cores=global$mc.cores.internal)
}

n.tot.ind <- function(ind, global) {
	if (is.null(ind$n.tot))
		sum(c(ind$gam1, ind$gam2))
	else
		ind$n.tot
}

n.tot.pop <- function(pop, global) {
	mcsapply(pop, n.tot.ind, global=global, mc.cores=global$mc.cores.internal)
}

n.piRNA.ind <- function(ind, global) {
	if (is.null(ind$n.piRNA))
		sum(c(which(ind$gam1), which(ind$gam2)) %in% global$piRNA.loci)
	else
		ind$n.piRNA
}

n.piRNA.pop <- function(pop, global) {
	mcsapply(pop, n.piRNA.ind, global=global, mc.cores=global$mc.cores.internal)
}


make.gamete <- function(ind) {
	ans <- ind$gam1
#~ 	rec <- sample(seq_along(ind$gam1), rbinom(1, length(ind$gam1), prob=0.5), replace=FALSE)
	rec <- runif(length(ind$gam1)) < 0.5
	ans[rec] <- ind$gam2[rec]
	ans
}

fitness.ind <- function(ind, global) {
	if (is.null(ind$fitness)) {
		n.sel <- n.sel.ind(ind, global)
		global$fitness.FUN(n.sel)
	} else 
		ind$fitness
}

fitness.pop <- function(pop, global) {
	sapply(pop, fitness.ind, global=global)
}

reproduction.ind <- function(par1, par2, global) {
	ind <- list(
		gam1=make.gamete(par1), 
		gam2=make.gamete(par2))
	class(ind) <- "ind"
	ind$n.tot   <- n.tot.ind(ind, global)
	ind$n.sel   <- n.sel.ind(ind, global)
	ind$n.piRNA <- n.piRNA.ind(ind, global)
	ind$fitness <- fitness.ind(ind, global)
	ind
}

reproduction.pop <- function(pop, global) {
	fitnesses <- fitness.pop(pop)
	N <- length(pop)
	parent1 <- sample(1:N, N, replace=TRUE, prob=fitnesses)
	parent2 <- sample(1:N, N, replace=TRUE, prob=fitnesses)
	newpop <- mclapply(1:N, function(i) { 
			reproduction.ind(pop[[parent1[i]]], pop[[parent2[i]]], global)
		}, mc.cores=global$mc.cores.internal)
	newpop
}

transposition.ind <- function(ind, global) {
	n.piRNA <- n.piRNA.ind(ind, global)
	n.tot <- n.tot.ind(ind, global)
	transp.rate <- global$u * global$regulation.FUN(n.piRNA)
	new.insertions <- rpois(1, n.tot*transp.rate)
	if (new.insertions > 0)
		for (nn in 1:new.insertions)
			if (runif(1) < 0.5)
				ind$gam1[sample(seq_along(ind$gam1), 1)] <- TRUE
			else 
				ind$gam1[sample(seq_along(ind$gam1), 1)] <- TRUE

	ind$n.tot   <- n.tot.ind(ind, global)
	ind$n.piRNA <- n.piRNA.ind(ind, global)
	ind$n.sel   <- n.sel.ind(ind, global)
	ind
}

transposition.pop <- function(pop, global) {
	mclapply(pop, transposition.ind, global=global, mc.cores=global$mc.cores.internal)
}


init.ind <- function(global) {
	ind <- list()
	class(ind) <- "ind"
	ind$gam1 <- ind$gam2 <- rep(FALSE, global$nb.loci)
	ind$gam1[sample(seq_along(ind$gam1), min(global$nb.loci, rpois(1, global$init.TE.ind/2)), replace=FALSE)] <- TRUE
	ind$gam2[sample(seq_along(ind$gam2), min(global$nb.loci, rpois(1, global$init.TE.ind/2)), replace=FALSE)] <- TRUE
	ind$n.tot <- n.tot.ind(ind, global)
	ind$n.sel <- n.sel.ind(ind, global)
	ind$n.piRNA <- n.piRNA.ind(ind, global)
	ind$fitness <- fitness.ind(ind, global)
	ind
}

init.pop <- function(global) {
	mclapply(1:global$N, function(i) init.ind(global), mc.cores=global$mc.cores.internal)
}

summary.pop <- function(pop, global) {
	ans <- list()
	n.tot <- n.tot.pop(pop, global)
	n.sel <- n.sel.pop(pop, global)
	n.piRNA <- n.piRNA.pop(pop, global)
	fitnesses <- fitness.pop(pop, global)
	ans$n.tot.mean <- mean(n.tot)
	ans$n.tot.var  <- var (n.tot)
	ans$n.sel.mean <- mean(n.sel)
	ans$n.sel.var  <- var (n.sel)
	ans$n.piRNA.mean<-mean(n.piRNA)
	ans$n.piRNA.var<- var (n.piRNA)
	ans$fitness.mean <- mean(fitnesses)
	ans$fitness.var  <- var (fitnesses)
	
	ans
}

run.one.simul <- function(global) {
	pop <- init.pop(global)
	ans <- data.frame(summary.pop(pop, global))
	for (gg in 1:global$G) {
		pop <- reproduction.pop(pop, global)
		pop <- transposition.pop(pop, global)
		ans <- rbind(ans, summary.pop(pop, global))
	}
	rownames(ans) <- as.character(0:global$G)
	ans
}

run.simul <- function(user.global=list()) {
	global <- global.default
	valid.user.global <- names(user.global)[names(user.global) %in% names(global)]
	if (length(valid.user.global) != length(user.global)) {
		warning(names(user.global)[! names(user.global) %in% global], ": unknown parameter names.")
	}
	global[valid.user.global] <- user.global[valid.user.global]
	if (is.na(global$mc.cores.internal))
		global$mc.cores.internal <- max(1, global$mc.cores %/% global$rep)
	mclapply(1:global$rep, function(i) run.one.simul(global=global), mc.cores=global$mc.cores)
}
