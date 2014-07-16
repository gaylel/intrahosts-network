args <- commandArgs(TRUE) ;
fname <- args[1]
load(fname)
library("coda")
library("ape")
dyn.load("intrahost_R2.so") 
source("ih_model2.R")
postscript("traj.eps")

t <- 3
bet = 2.7e-5
p=1.2e-2
c=3.0
delt=4.0
V0 = 9.3e-2
T0 = 4e8
mr=5.5e-6
nsites = 903
ns = NULL

params <- list(bet= bet, p=p, c=c, delt=delt, V0=V0, T0=T0,  nsites=nsites, Ns = ns, t=t, ss= 0.01, mr=mr)
traj <- ih_model_plot3(mcmc.out, params)
traj[is.na(traj)] <- 0
ti <- traj[nrow(traj), ]
traj <- (traj[-nrow(traj), ])
traj[is.na(traj)] <- 0

q <- apply((traj), 2, quantile, c(0.25, 0.5, 0.75), na.rm=TRUE)
plot(ti[500:1000], q[3, 500:1000], "l", xlab="t", ylab="i(t)", lty=3, lwd=3)
lines(ti[500:1000], q[2, 500:1000], "l", lwd=3)
lines(ti[500:1000], q[1, 500:1000], "l", lty=3, lwd=3)
#plot(exp(mcmc.out))
dev.off()
library("OutbreakTools")
source("ih_io.R")
source("ih_model2.R")


x <- ih_io_readdata() 

# calculate unique sequences
nsites <- length(get.dna(x)[[1]][1,])
uniqseq <- dna2uniqSequences(get.dna(x)[[1]])
IDcounts <- do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts <- as.data.frame(IDcounts[order(-IDcounts[, 1]), ])
colnames(IDcounts) <- c("count")

# consensus sequence
cseq <- row.names(IDcounts)[1]

inds <- get.individuals(x)
inds <- c(1)
sfs <- ih_io_sfs(x, uniqseq, cseq, inds)
postscript("sfs.eps")
barplot(sfs[1:20])
dev.off()
