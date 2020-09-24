

cc83 <- function(u=function(n) 0.1, v=0, n0=1, Tmax=100, dlw=function(n) -0.01*n) {
	ans <- c(n0, rep(NA, Tmax))
	for (t in 1:Tmax) {
		ans[t+1] <- ans[t] + ans[t]*(u(ans[t]) - v + ans[t]*dlw(ans[t]))
	}
	ans
}

mod1.nosel <- function(u, rho, v=0, n0=1, p0=0, Tmax=100) {
	ans <- rbind(c(n0, rep(NA, Tmax)), c(p0, rep(NA, Tmax)))
	rownames(ans) <- c("n","p")
	for (t in 1:Tmax) {
		ans["n",t+1] <- ans["n",t] + ans["n",t]*(u*(1-ans["p",t])^2 - v)
		ans["p",t+1] <- ans["p",t] + u*rho*ans["n",t]*(1-ans["p",t])^2
	}
	ans
}

mod1.timetoeq <- function(..., Tmax=10000, thresh=0.99) {
	aa <- do.call(mod1.nosel, c(list(Tmax=Tmax), list(...)))
	list(t=which(aa["p",] > thresh)[1], n=aa["n",Tmax])
}
