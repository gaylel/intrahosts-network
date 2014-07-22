args <- commandArgs(TRUE) ;
inds <- args[1]
paramfile <- args[2]
outdir <- args[3]

#inds <- 26
#paramfile <- "~/Documents/git/intrahosts-network/run/20140721/opt.params"  
#outdir <- "~/Documents/git/intrahosts-network/run/20140721/"

library("OutbreakTools")
library("coda")

#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model3.R")
source(paramfile)
datafile <- paste(outdir, "/intrahost.", inds, ".RData", sep="")
if (file.exists(datafile)
{
load(datafile)
trajfile <- paste(outdir, "/mcmc.", inds, ".pdf", sep="")
esfile <- paste(outdir, "/mcmc.", inds, ".esize", sep="")
par_init <- ih_model_getparams(params.isopt, params.init)
#traj <- ih_model_plot3(mcmc.out, params.init, par_init)

params <- c(params.init)
print(params)

traj <- ih_model_plottraj(mcmc.out, params, names(par_init))

traj[is.na(traj)] <- 0
ti <- traj[nrow(traj), ]
traj <- (traj[-nrow(traj), ])
traj[is.na(traj)] <- 0

q <- apply((traj), 2, quantile, c(0.25, 0.5, 0.75), na.rm=TRUE)

pdf(trajfile)
colnames(mcmc.out) <- names(par_init)
plot(exp(mcmc.out))

plot(ti[500:1000], q[3, 500:1000], "l", xlab="t", ylab="i(t)", lty=3, lwd=3)
lines(ti[500:1000], q[2, 500:1000], "l", lwd=3)
lines(ti[500:1000], q[1, 500:1000], "l", lty=3, lwd=3)
lines(D$ts - D$ts[1], D$cn * params$alph, "p")
#plot(exp(mcmc.out))
dev.off()

write.table(effectiveSize(exp(mcmc.out)), file=esfile, quote=FALSE, col.names=FALSE)
}