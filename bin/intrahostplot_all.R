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
cn <- NULL
for (inds in seq(1, 50))
{
	datafile <- paste(outdir, "/intrahost.", inds, ".RData", sep="")
	if (file.exists(datafile) && inds!=27)
	{
		load(datafile)
		mcmc.out <- exp(mcmc.out)
		#med <- rbind(med, apply((mcmc.out), 2, median))
		med <- rbind(med, mcmc.out)
		if (length(D$cn) == 1)		
		cn <- rbind(cn, c(apply(mcmc.out , 2, median), D$cn, D$S, D$Ns))
	}
}

colnames(med) <- names(par_init)
histfile <- paste(outdir, "/mcmc.hist.pdf", sep="")

pdf(histfile)
for (i in seq(1, ncol(med)))
hist(med[,i], main=colnames(med)[i], xlab="", 100)
#print(med)

colnames(cn) <- c(names(par_init), "cn", "S", "Ns")

plot(cn[,"cn"], cn[,"t"])
plot(cn[,"S"]/cn[,"Ns"], cn[,"t"])
plot(cn[,"S"], cn[,"V0"])
print(cn)
dev.off()

