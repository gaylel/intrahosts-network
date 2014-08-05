args <- commandArgs(TRUE) ;
testno <- args[1]
expdir <- args[2]
outdir <- args[3]

#inds <- 26
#paramfile <- "~/Documents/git/intrahosts-network/run/20140721/opt.params"  
#outdir <- "~/Documents/git/intrahosts-network/run/20140721/"

library("OutbreakTools")
library("coda")
library("ape")
#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model3.R")

			
			
sumofexp <- function(loga, sort=0)
{
	## returns log(sum(exp(loga)))
	if (sort==1)
		loga<-sort(loga)
	
	L <- loga[1]	
	if (length(loga) > 1)
	{
	for (i in seq(2, length(loga)))
	{
		Lv <- c(L,loga[i])
		mll <- max(Lv)
		L <- log(sum(exp(Lv - mll))) + mll 
	}
	}
	#L <- L- log(length(loga))
	return(L)
}



expparamfile <- paste(expdir, "/exp.params", sep="")
paramfile <- paste(expdir, "/opt.params", sep="")
source(paramfile)
source(expparamfile)
allvars <- expand.grid(c(exp.params$dvars, exp.params$pvars))
varnames <- names(allvars)
	
switch(as.character(testno), 
	"1"={
		t_i <- seq(2,6,by=0.1)
		#for (v_it in seq(1, nrow(allvars)))
		#{
			source(paramfile)

		#	# output files
		#	varstr <- NULL
		#	for (i in seq(1, length(varnames)))
		#	{
		#		varstr <- paste(varstr, paste(varnames[i], "_", allvars[v_it, i], "_", sep=""), sep="")	
		#	}

		#	outdir <- paste(expdir, "/", varstr, sep="")
					

	
	
	
			llpdffile <- paste(outdir, "/ll_t_post.pdf", sep="")
			ll_all <- NULL
			llv <- read.table(paste(outdir, "/ll_v_t.dat", sep=""))
			llv <- llv[,2:ncol(llv)]
			
			
			for (i in exp.params$drep) 
			{
				L <- sumofexp(unlist(llv[i,])) ; 
				ll_all <- rbind(ll_all, exp(unlist(llv[i,] - L ))) 
			}
			
			pdf(llpdffile)
			plot(t_i, q[2,], type="l", lwd = 2, xlab="t", ylab="p(t|D)", ylim=c(0,0.5))
			#polygon(c(t_i, rev(t_i)), c(q[1,], rev(q[3,])), col="gray", border="gray")
			#lines(t_i, q[2,], type="l", lwd = 2)
			dev.off()
			
			
	}
)