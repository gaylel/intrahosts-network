# fitting of single intrahost model based on site frequency spectrum
mcmc.ll <- NULL

ih_get_traj_R <- function(bet, p, c, delt, v0, T0, ss)
{
	traj <- .Call("ih_get_traj", bet, p, c, delt, v0, T0, ss)
	return(traj)
}

ih_model_init <- function(t, Ns)
{
	# time point	t
	# number of sequences Ns
	
}

ih_model_fit <- function(params)
{
	bet <- params$bet
	c <- params$c
	p <- params$p
	delt <- params$delt
	V0 <- params$V0
	T0 <- params$T0
	t <- params$t
	Ns <- params$Ns
	ss <- 0.01
	ts <- c(10, 12)
	#ih <- .Call("ih_R_drawtraj", bet, p, c, delt, V0, T0, ss, Ns, t, 1)
	ih <- .Call("ih_R_drawtraj2", bet, p, c, delt, V0, T0, ss, Ns, ts, t, 1)
	return(ih)
}

ih_model_sfsll <- function(params, D)
{
	bet <- params$bet
	c <- params$c
	p <- params$p
	delt <- params$delt
	V0 <- params$V0
	T0 <- params$T0
	t <- params$t
	Ns <- params$Ns
	ss <- 0.01
	alph = 1000
	ih <-.Call("ih_R_getsfs", bet, p, c, delt, V0, T0, ss, t, D)
	#print(ih)
	ll <- log(0)
	if (ih$err != 0)
	{
		ll <- ih$ll
		ll3 <- 0
		r <- params$nsites * ih$ES * params$mr
		ll2 <- dpois(D$S, r, log=TRUE)
		#cn = c(1090, 2310)
		#cn = c(1370, 75700)
		cn = D$cn
		for (i in seq(1, length(ih$cn)))
		{
		#ll3 <- ll3 + dnorm(log10(10000 * cn[i]), log10(ih$cn[i]), 1, log = TRUE)
		}
		ll <- ll + ll2 + ll3
	
	# add priors
	
	}
	return(ll) 
}



ih_model_opt <- function(pars, params, D)
{
	t <- exp(pars[1])
	bet <- exp(pars[2])
	c <- exp(pars[3])
	p <- exp(pars[4])
	delt <- exp(pars[5])
	V0 <- exp(pars[6])
	#mr <-exp(pars[6])
	
	params$t <- t
	params$bet <- bet
	params$c <- c
	params$delt <- delt
	params$p <- p
	params$V0 <- V0
	#params$mr <- mr
	#ih <- ih_model_fit(params)
	ll <- log(0)
	ll <- ih_model_sfsll(params, D)
	#if (ih$ET[1] >= 0)
	#{
	#ll <- ih_model_ll(ih, data$sfs, data$S, params)
	# add priors
	#ll<- ll + dgamma(pars, shape=1e-5, scale=1e5,log=TRUE)
	#}
	ll<- ll + dgamma(pars, shape=1e-5, scale=1e5,log=TRUE)
	ll <- ll + dgamma(t, shape=40, scale=0.1, log=TRUE)
	#mcmc.ll <<- c(mcmc.ll, ll)
	return(ll) 
}

ih_ancestor_function <- function(ih, t)
{
	nseq <- length(ih$ET) + 1
	bt <- rev(t - cumsum(rev(ih$ET)))
	At <- c(0, bt, t)
	At <- rbind(c(seq(1,nseq),nseq), At)
	return(At)
}

ih_model_ll <- function(ih, counts, S, params)
{
	ll <- log(0)
	ll1 <- ll
	ll2 <- ll
	#print(ih)
	if (ih$err == 1)
	{
		ll1 <- dmultinom(counts, prob=ih$ESS, log = TRUE)
	
	r <- params$nsites * ih$ES * params$mr
	ll2 <- dpois(S, r, log=TRUE)
	ll3 <- dnorm(log10(2720000), log10(ih$cn), 1, log = TRUE)
	ll <- ll1 + ll2 + ll3
	}
	#print(ih$ES)
	#print(sprintf("ll=%8.4f, mu=%8.4f, ll(mu)=%8.4f, ll1=%8.4f, t=%8.4f", ll, nsites * params$mr *ih$ES, ll2, ll1, params$t ))
	#print(ih$cn)
	return(ll)
}

ih_model_mcmc <- function(par_init, params, mcmc.params, D)
{
	# optimise parameters in logspace.
	par_init <- log(par_init)
	#V[2, 3] <- 0.5
	#V[3, 2] <- 0.5
	#mcmc.out <- MCMCmetrop1R(ih_model_opt, theta.init = par_init, params=params, data=data, mcmc.params)
	mcmc.out <- do.call(MCMCmetrop1R, c(list(fun=ih_model_opt, theta.init = par_init, params=params, D=D), mcmc.params))
	return(list(mcmc.out=mcmc.out, mcmc.ll=mcmc.ll))
}

ih_model_gridsearch <- function(params, data)
{
	t_in <- seq(7, 14, by=0.1)
	gam_in <- seq(0.1, 0.5, by=0.02)
	bet_in <- seq(1,3,by=0.02)
	ll_all <- array(-log(0), dim=c(length(t_in), length(bet_in), length(gam_in)))

	for (i in seq(1, length(t_in)))
	{
		for (j in seq(1, length(bet_in)))
		{
			for (k in seq(1, length(gam_in)))
			{
				if (bet_in[j] >= gam_in[k])
					ll_all[i,j,k] <- ih_model_opt(c(t_in[i], bet_in[j], gam_in[k]), params, data)
			}
		}
	}
	
	

	save(ll_all, gam_in, t_in, bet_in, file="ll_all.RData")
}



ih_model_plot <- function(ih, t)
{
	postscript(file="ih_plot.eps", onefile=TRUE, horizontal=TRUE)
	par(mfrow=c(2,2))
	plot(ih$sir$t, ih$sir$i, "l")	
	barplot(ih$ESS)
	At <- ih_ancestor_function(ih, t)
	plot(At[2,], At[1,], "s", xlim=c(0, max(ih$sir$t)))
	dev.off()
}

ih_model_plot2 <- function(mcmc.in, params)
{
	Tmax = 20
	N <- 1000
	L <- nrow(mcmc.in)
	ti <- seq(-Tmax, 0, length.out=N+1)
	traj <- NULL
	print(L)
	for (i in seq(1, L, by=100))
	{
		params$t <- exp(mcmc.in[i, 1])
		params$bet <- exp(mcmc.in[i, 2])
		params$c <- exp(mcmc.in[i, 3])
		params$p <- exp(mcmc.in[i, 4])
		params$delt <- exp(mcmc.in[i, 5])
		ih <- ih_model_fit(params)
		T <- ih$sir$t - params$t
		ap <- approx(T, ih$sir$i, ti, method="constant")
		traj <- rbind(traj, ap$y)
	}
	traj <- rbind(traj, ti)
	return(traj)
}


ih_model_plot3 <- function(mcmc.in, params)
{
	Tmax = 20
	N <- 1000
	L <- nrow(mcmc.in)
	ti <- seq(-Tmax, 0, length.out=N+1)
	traj <- NULL
	print(mcmc.in)
	for (i in seq(1, L, by=10))
	{
		params$t <- exp(mcmc.in[i, 1])
		params$bet <- exp(mcmc.in[i, 2])
		params$c <- exp(mcmc.in[i, 3])
		params$p <- exp(mcmc.in[i, 4])
		params$delt <- exp(mcmc.in[i, 5])
		params$V0 <- exp(mcmc.in[i, 6])
		sir <- ih_get_traj_R(params$bet, params$p, params$c, params$delt, params$V0, params$T0, params$ss)
		sir <- sir[[1]]
		T <- sir$t - (params$t)
		ap <- approx(T, (sir$v), ti, method="constant")
		traj <- rbind(traj, ap$y)
	}
	traj <- rbind(traj, ti)
	return(traj)
}

