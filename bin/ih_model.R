# fitting of single intrahost model based on site frequency spectrum

ih_model_init <- function(t, Ns)
{
	# time point	t
	# number of sequences Ns
	
}

ih_model_fit <- function(t, Ns, params)
{
	bet <- params$bet
	gam <- params$gam
	I0 <- params$I0
	NS <- params$NS
	
	ss <- 0.01
	
	ih <- .Call("ih_R_drawSIR", bet, gam, I0, NS, ss, Ns, t)
	return(ih)
}

ih_ancestor_function <- function(ih, t)
{
	nseq <- length(ih$ET) + 1
	bt <- rev(t - cumsum(rev(ih$ET)))
	At <- c(0, bt, t)
	At <- rbind(c(seq(1,nseq),nseq), At)
	return(At)
}

ih_model_ll <- function(ih, counts)
{
	ll <- dmultinom(counts, prob=ih$ESS, log = TRUE)
	return(ll)
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