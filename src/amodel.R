source("../src/sim_functions.R")

simulator.default = "Rsimulator"
simulator.default = "Simulicron"

cc83 <- function(u=function(n) 0.1, v=0, n0=1, Tmax=100, dlw=function(n) -0.01*n) {
	ans <- c(n0, rep(NA, Tmax))
	for (t in 1:Tmax) {
		ans[t+1] <- ans[t] + ans[t]*(u(ans[t]) - v + ans[t]*dlw(ans[t]))
	}
	ans
}

model <- function(t, y, parms) {
	with(as.list(c(y, parms)), {
		dn <- n*u*(1-p)^(2*k) - n*s*(1+2*u)
		dp <- n*u*pi/k*(1-p)^(2*k) - sp*p*(1-p)/(1-sp*p)
		list(c(dn, dp))
	})
}

amodel <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, Tmax=100) {
	library(deSolve)
		
	init <- setNames(c(n0, p0), nm=c("n","p"))
	params <- setNames(c(u=u, pi=pi, k=k, s=s, sp=sp), nm=c("u","pi","k","s","sp"))
	times <- seq(0, Tmax, 1)
	
	oo <- ode(y=init, times=times, func=model, parms=params)

	ans <- list(
		n=oo[,"n"],
		p=oo[,"p"])
	ans
}


pred.eq <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, r=NA) {
	if (p0 != 0) warning("Most models assume that p0=0.")
	ss <- s*(1+2*u)
	pp <- 1-(ss/u)^(1/2/k)
	if (pp < 0) pp <- 0
	if (s==0) {
		return(list(Eq=list(
						n=n0+k/pi,
						p=1)))
	} else { #s != 0
		if (sp == 0) {
			nn <- n0 + (k/pi)*(1- (2*k*(ss/u)^(1/2/k) - ss/u)/(2*k-1))
			return(list(Max=list(
							n = nn,
							p = pp),
						Eq=list(
							n = 0, 
#~ 							p = 1 - ((1-pp)^(1-2*k) + (2*k-1)*u*pi*nn/k/ss)^(1/(1-2*k)) )
							p = 1- ((u/ss)*(2*k-1)*pp +1) ^ (1/(1-2*k)) )
#~ 							p = pp - sqrt(abs((1/k)*(ss/u)^(1/2/k)*(2/(1-2*k)*(ss/u)^(1/2/k) - pi*nn/k))) )
						))
			} else { #sp != 0
				return(list(Eq=list(
#~ 									n=k*sp/(pi*u)*((s/u)^(1/2/k)+(s/u)^(-1/2/k) - 2),
#~ 									n=k*s/(pi*u)*(1-(s/u)^(1/2/k))/((s/u)^((2*k-1)/2/k)*(1-s*(s/u)^(1/2/k))),
									n=k*s*pp*(ss/u)^(1/2/k)/pi/ss/(1-s*pp),
									p=pp)))
			}
	}
}

num.eq <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0) {
	library(deSolve)
	library(phaseR)
		
	init <- setNames(c(n0, p0), nm=c("n","p"))
	params <- setNames(c(u=u, pi=pi, k=k, s=s, sp=sp), nm=c("u","pi","k","s","sp"))
	times <- seq(0, 1000, 1)

	oo <- ode(y=init, times=times, func=model, parms=params)
	eq <- findEquilibrium(model, y0=oo[nrow(oo),-1], parameters=params, state.names=names(init), summary=FALSE)

	return(list(Eq=list(n=eq$ystar["n",1], p=eq$ystar["p",1])))
}

jacob.dtdc <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, ...) {
	ss <- s*(1+2*u)
	pp <- 1-(ss/u)^(1/2/k)
	nn <- k*s*pp*(ss/u)^(1/2/k)/pi/ss/(1-s*pp)
	rbind(
		  c(0,       -2*k*nn*u*(ss/u)^((2*k-1)/2/k)),
		  c(pi*ss/k, -2*nn*pi*u*(ss/u)^((2*k-1)/2/k) + (1-s)/(1-s*pp)^2 - 1))
}

simmodel <- function(u=0.1, pi=0.03, s=0, k=1, sp=0, n0=1, p0=0, r=0, N=10000, Tmax=100, rep=1, use.cache=TRUE, cache.dir="../cache/", mean=TRUE, simulator=simulator.default) {
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
		init.TE.ind       = n0,
		simulator         = simulator
	)
	simpar$regulation.FUN <- function(n.piRNA) if (n.piRNA == 0) 1 else 0
	# This tries to go around R lazy evaluation of closures ... what a mess
	ffs <- paste0('function(n.sel) exp(-' ,eval(s), '*n.sel)')
	simpar$fitness.FUN <-  eval(parse(text=ffs))
	
	ffr <- paste0('function(n.piRNA, n) (if (n.piRNA == 0) 1 else 0)/(1+', eval(r), '*n)')
	simpar$regulation.FUN <- eval(parse(text=ffr))
	
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
					file.path(cache.dir, paste0(paste(digest(capture.output(dput(simpar))), collapse=""), ".rds")) 
				else 
					NULL
	if (use.cache && file.exists(file.name)) {
		simres <- readRDS(file.name)
	} else {
		simres <- run.simul(simpar, simulator=simulator)
	}
	if (use.cache && !file.exists(file.name)) {
		saveRDS(simres, file.name)
	}
	if (mean) {
		ans <- list(
			n=rowMeans(sapply(simres, function(i) i$n.tot.mean)),
			p=rowMeans(sapply(simres, function(i) i$n.piRNA.mean))/2/k)
	} else {
		ans <- list(
			n=sapply(simres, function(i) i$n.tot.mean),
			p=sapply(simres, function(i) i$n.piRNA.mean)/2/k)
	}
	names(ans$n) <- names(ans$p) <- rownames(simres[[1]])
	ans
}

plot.model.dyn <- function(model.default, model.par, what="n", pred=TRUE, sim=FALSE, legend=TRUE, Tmax=100, N=10000, rep=1, nb.simpt=21, use.cache=TRUE,
	xlab="Generations", ylab=if(what=="n") "Copy number" else "Cluster frequency", xlim=c(0,Tmax), ylim=NA,
	legend.pos="topleft", simulator=simulator.default, ...) {

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
			do.call(simmodel, c(as.list(pp), list(N=N, Tmax=Tmax, rep=rep, use.cache=use.cache, simulator=simulator)))
		})
	}
	
	if (is.na(ylim)) ylim <- c(0, 1.2*max(unlist(sapply(dyn.res, "[[", what)), unlist(sapply(pred.res, "[[", what)), unlist(sapply(sim.res, "[[", what))))
	
	plot(NULL, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)

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
