args <- commandArgs(TRUE) ;
inds <- args[1]
paramfile <- args[2]
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
source(paramfile)
datafile <- paste(outdir, "/intrahost.", inds, ".RData", sep="")
if (file.exists(datafile))
{
load(datafile)
trajfile <- paste(outdir, "/mcmc.", inds, ".pdf", sep="")
esfile <- paste(outdir, "/mcmc.", inds, ".esize", sep="")
btfile <- paste(outdir, "/mcmc.", inds, ".ebt", sep="")
par_init <- ih_model_getparams(params.isopt, params.init)
#traj <- ih_model_plot3(mcmc.out, params.init, par_init)

params <- c(params.init)
print(params)

traj <- ih_model_plottraj(mcmc.out, params, names(par_init), -7, 5)

traj[is.na(traj)] <- 0
ti <- traj[nrow(traj), ]
traj <- (traj[-nrow(traj), ])
#traj[is.na(traj)] <- 0
#plot(traj[1,])
q <- apply((traj), 2, quantile, c(0.25, 0.5, 0.75), na.rm=TRUE)

pdf(trajfile)
colnames(mcmc.out) <- names(par_init)
plot(exp(mcmc.out))

plot(ti, q[3, ], "l", xlab="t", ylab="v(t)", lty=3, lwd=3)
lines(ti, q[2, ], "l", lwd=3)
lines(ti, q[1, ], "l", lty=3, lwd=3)
lines(D$ts - D$ts[1], log10(D$cn * params$alph), "p")
#plot(exp(mcmc.out))

write.table(effectiveSize(exp(mcmc.out)), file=esfile, quote=FALSE, col.names=FALSE)



ebt <- ih_model_get_ebt_all(mcmc.out, params, names(par_init), D)
write.table(ebt, file=btfile, quote=FALSE, col.names=FALSE)
#print(ebt)
#q <- apply(ebt, 2, quantile, c(0.25, 0.5, 0.75), na.rm=TRUE)



boxplot(ebt, outline=FALSE)
tr2 <- read.tree(paste(outdir, "/tr.nwk", sep=""))

tr2_c <- coalescent.intervals(tr2)
lines(tr2_c$lineages - 1, tr2_c$interval.length, "l", col="red", lwd=2)
dev.off()
}
