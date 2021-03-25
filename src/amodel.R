source("../src/sim_functions.R")

cc83 <- function(u=function(n) 0.1, v=0, n0=1, Tmax=100, dlw=function(n) -0.01*n) {
	ans <- c(n0, rep(NA, Tmax))
	for (t in 1:Tmax) {
		ans[t+1] <- ans[t] + ans[t]*(u(ans[t]) - v + ans[t]*dlw(ans[t]))
	}
	ans
}

amodel <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, Tmax=100, corr=FALSE) {
	ans <- list(
		n=c(n0, rep(NA, Tmax)),
		p=c(p0, rep(NA, Tmax)))
	for (t in 1:Tmax) {
		
#~ 		# selection
#~ 		p <- ans$p[t]
#~ 		n <- ans$n[t]
#~ 		ans$n[t+1] <- n - s*n*(1+2*u)
#~ 		ans$p[t+1] <- p - sp*p*(1-p)/(1-sp*p)
		
#~ 		# transposition
#~ 		p <- ans$p[t+1]
#~ 		n <- ans$n[t+1]
#~ 		cf <- (1-p)^(2*k)		
#~ 		ans$n[t+1] <- n + n*u*cf
#~ 		ans$p[t+1] <- p + n*u*pi*cf/k 

		p <- ans$p[t]
		n <- ans$n[t]
		cf <- (1-p)^(2*k)		
		ans$n[t+1] <- n + n*u*cf - n*s*(1+2*u)
		ans$p[t+1] <- p + n*u*pi*cf/k - sp*p*(1-p)/(1-sp*p)
	}
	ans
}


pred.eq <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, corr=TRUE) {
	if (corr)
		s <- s*(1+4*u)
	if (p0 != 0) warning("Most models assume that p0=0.")
	
	if (s==0) {
		return(list(Eq=list(
						n=n0+k/pi,
						p=1)))
	} else { #s != 0
		if (sp == 0) {
			return(list(Max=list(
							n= n0 + (k/pi)*(1-(s/u)^(1/2/k) + s/u/(2*k-1)*(1-(s/u)^(1/2/k-1))),
#~ 							n=n0+(k/pi)*(1-(s/u)^(1/2/k))^2,
							p=1-(s/u)^(1/(2*k))),
						Eq=list(
							n=0,
							p=1-1/(u*(1+n0*pi)/s-1))))
			} else { #sp != 0
				return(list(Eq=list(
#~ 									n=k*sp/(pi*u)*((s/u)^(1/2/k)+(s/u)^(-1/2/k) - 2),
#~ 									n=k*s/(pi*u)*(1-(s/u)^(1/2/k))/((s/u)^((2*k-1)/2/k)*(1-s*(s/u)^(1/2/k))),
									n=(k/pi)*((s/u)^(1/k)-(s/u)^(1/2/k))/(s*(s/u)^(1/2/k)-1),
									p=(1-(s/u)^(1/2/k)))))
			}
	}
}

simmodel <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, N=10000, Tmax=100, rep=1, use.cache=TRUE, cache.dir="../cache/") {
	simpar <- list(
		nb.loci           = 1000,
		neutral.loci      = if (sp==0) 1:k else numeric(0), 
		piRNA.loci        = 1:k,
		piRNA.prob        = pi,
		u                 = u,
		G                 = Tmax,
		N                 = N,
		rep               = rep,
		summary.every     = 1,
		init.TE.ind       = n0
	)
	simpar$regulation.FUN <- function(n.piRNA) if (n.piRNA == 0) 1 else 0
	# This tries to go around R lazy evaluation of closures ... what a mess
	ffs <- paste0('function(n.sel) exp(-' ,eval(s), '*n.sel)')
	simpar$fitness.FUN <-  eval(parse(text=ffs))
	
	if (use.cache && is.null(cache.dir)) {
		warning("Impossible to use the cache (no cache directory provided).")
		use.cache <- FALSE
	}
	if (use.cache && !require(digest, quietly=TRUE)) {
		warning("Impossible to use the cache (library digest not available).")
		use.cache <- FALSE
	}
	if (use.cache && !dir.exists(cache.dir)) {
		# This would deserve a better path management
		dir.create(cache.dir) 
	}
	file.name <- if (use.cache) 
					file.path(cache.dir, paste0(digest(capture.output(dput(simpar))), ".rds")) 
				else 
					NULL
	
	if (use.cache && file.exists(file.name)) {
		simres <- readRDS(file.name)
	} else {
		simres <- run.simul(simpar)
	}
	if (use.cache) {
		saveRDS(simres, file.name)
	}
	ans <- list(
		n=rowMeans(sapply(simres, function(i) i$n.tot.mean)),
		p=rowMeans(sapply(simres, function(i) i$n.piRNA.mean))/2/k)
	names(ans$n) <- names(ans$p) <- rownames(simres[[1]])
	ans
}

plot.model.dyn <- function(model.default, model.par, what="n", pred=TRUE, sim=FALSE, legend=TRUE, Tmax=100, N=10000, rep=1, nb.simpt=21, use.cache=TRUE,
	xlab="Generations", ylab=if(what=="n") "Copy number" else "Cluster frequency", xlim=c(0,Tmax), ylim=NA,
	legend.pos="topleft") {

	dyn.res <- lapply(model.par, function(mm) { 
		pp <- model.default
		pp[names(mm)] <- mm
		do.call(amodel, c(as.list(pp), list(Tmax=Tmax)))
	})
	
	if (pred){
		pred.res <- lapply(model.par, function(mm) {
			pp <- model.default
			pp[names(mm)] <- mm				
			do.call(pred.eq, as.list(pp))
		})
	}
	
	if (sim) {
		sim.res <- lapply(model.par, function(mm) {
			pp <- model.default
			pp[names(mm)] <- mm				
			do.call(simmodel, c(as.list(pp), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache )))
		})
	}
	
	if (is.na(ylim)) ylim <- c(0, 1.2*max(unlist(sapply(dyn.res, "[[", what)), unlist(sapply(pred.res, "[[", what)), unlist(sapply(sim.res, "[[", what))))
	
	plot(NULL, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)

	for (ip in seq_along(model.par)) {
		lines(0:Tmax, dyn.res[[ip]][[what]], col=col[ip])
		if (pred) {
			if ("Max" %in% names(pred.res[[ip]]))
				abline(h=pred.res[[ip]]$Max[[what]], lty=lty.max, col=col[ip])
			if ("Eq" %in% names(pred.res[[ip]]))
				abline(h=pred.res[[ip]]$Eq[[what]], lty=lty.eq, col=col[ip])	
		}
		if (sim) {
			xx <- as.numeric(names(sim.res[[ip]][[what]]))
			xpl <- unique(seq(1, length(xx), length.out=nb.simpt))
			points(xx[xpl], sim.res[[ip]][[what]][xpl], pch=1, col=col[ip])
		}
	}
	
	if (legend) {
		legend.labels <- sapply(seq_along(model.par), function(ip) {paste0(names(model.par[[ip]]), "=", model.par[[ip]], collapse=", ")})
		legend(legend.pos, lty=1, col=col[seq_along(model.par)], legend=legend.labels, bty="n")
	}
}
