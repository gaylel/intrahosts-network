args <- commandArgs(TRUE) ;
inds <- args[1]

library("OutbreakTools")
library("ape")
library("MCMCpack")

#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model2.R")


#inds <- get.individuals(x)
#inds <- c(25)
tstart <- 20

####################################################################################
# Reading in horseflu data, calculcating (joint) sfs, and copy numbers

cn_fname <- "../data/journal.ppat.1003081.s005.csv"

# read in horseflu data
x <- ih_io_readdata() 
# read in copy number data
vcn <- ih_io_readcndata(cn_fname)

# calculate unique sequences
nsites <- length(get.dna(x)[[1]][1,])
uniqseq <- dna2uniqSequences(get.dna(x)[[1]])
IDcounts <- do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts <- as.data.frame(IDcounts[order(-IDcounts[, 1]), ])
colnames(IDcounts) <- c("count")



#sfs <- ih_io_sfs(x, uniqseq, cseq, inds)
#print(sfs)

seqs.meta <- subset(x, individual=inds)@dna@meta	
seqs.sids <- unique(seqs.meta$sampleID)
seqs.dates <- seqs.meta[match(seqs.sids, seqs.meta$sampleID), "date"]
o <- order(seqs.dates)
seqs.dates <- seqs.dates[o]
seqs.sids <- seqs.sids[o]

# get consensus sequence
cseq <- ih_get_consensus(x, uniqseq, seqs.sids)
# consensus sequence
#cseq <- row.names(IDcounts)[1]


# number of sequences input
ns <- numeric(length(seqs.sids)) 
for (i in seq(1, length(seqs.sids)))
{
	ns[i] <- length(which(seqs.meta[,"sampleID"]==seqs.sids[i]))
}

# dates input
ts <- ih_io_convertdates(seqs.dates, tstart)

# joint site frequency spectrum input
jsfs <- ih_io_jsfs(x, uniqseq, cseq, seqs.sids)

# copy number 
cn <- vcn$Copy.No[match(seqs.sids, vcn$Isolate)]

# input parameters
t <- 3
bet = 2.7e-5
p=1.2e-2
c=3.0
delt=4.0
V0 = 9.3e-2
T0 = 4e8
mr=5.5e-6

par_init <-c(t,bet,c,p,delt, V0) 
params <- list(bet= bet, p=p, c=c, delt=delt, V0=V0, T0=T0,  nsites=nsites, Ns = ns, t=t, mr=mr)
V <- diag(length(par_init))
mcmc.params <- list(tune=c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05), V=V, verbose=10, burnin=5000, mcmc=10000)
#data <- list(sfs=sfs[1:(ns[1]-1)], S=S)
D <- list(mut=jsfs$mut, count=jsfs$count, S=jsfs$S, ts=ts, Ns=ns, cn=cn)

#########################################################################################

# run model
mcmc.out <- ih_model_mcmc(par_init, params, mcmc.params, D)
ll <- mcmc.out$mcmc.ll
mcmc.out <- mcmc.out$mcmc.out
traj <- ih_model_plot2(mcmc.out, params)
#save(mcmc.out, traj, ll, file="mcmc_out_test.RData")


save(mcmc.out, traj, ll, file=paste("intrahost", inds, ".RData", 
sep=""))

