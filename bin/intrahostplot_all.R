args <- commandArgs(TRUE) ;
paramfile <- args[1]
outdir <- args[2]

library("OutbreakTools")
library("coda")

#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model3.R")
source(paramfile)
par_init <- ih_model_getparams(params.isopt, params.init)
med <- NULL
for (inds in seq(1, 50))
{
	datafile <- paste(outdir, "/intrahost.", inds, ".RData", sep="")
	if (file.exists(datafile))
	{
		load(datafile)
		colnames(mcmc.out) <- names(par_init)

		mcmc.out <- exp(mcmc.out)
		med <- rbind(med, apply((mcmc.out), 2, median))
	}
}

histfile <- paste(outdir, "/mcmc.hist.pdf", sep="")

pdf(histfile)
for (i in seq(1, ncol(mcmc.out)))
hist(mcmc.out[,i], main=colnames(mcmc.out)[i], xlab="", 50)
dev.off()
