args <- commandArgs(TRUE) ;
paramfile <- args[1]
indir <- args[2]
outdir <- args[3]
inds <- args[4]

library("OutbreakTools")
library("ape")
library("MCMCpack")

#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model3.R")


#inds <- get.individuals(x)
#inds <- c(25)
tstart <- 20

####################################################################################
# Reading in testdata, calculcating (joint) sfs, and copy numbers

x <- read.dna(paste(indir, "/seq.dat", sep=""))

#vcn <- ih_io_readcndata(cn_fname)

# calculate unique sequences
nsites <- length(x[1,])
uniqseq <- dna2uniqSequences(x)
IDcounts <- do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
IDcounts <- as.data.frame(IDcounts[order(-IDcounts[, 1]), ])
colnames(IDcounts) <- c("count")

# get consensus sequence
ids <- rownames(x)
cseq <- ih_consensus(uniqseq, ids)


# number of sequences input
ns <- length(ids)

# dates input
ts <- tstart

# joint site frequency spectrum input
jsfs <- ih_io_jsfs2(x, uniqseq, cseq)

# copy number 
cn <- 1000 #vcn$Copy.No[match(seqs.sids, vcn$Isolate)]

#### optimisation ########


# input parameters
source(paramfile)
par_init <- ih_model_getparams(params.isopt, params.init)
params <- c(params.init, list(nsites=nsites, Ns=ns))
#data <- list(sfs=sfs[1:(ns[1]-1)], S=S)
D <- list(mut=jsfs$mut, count=jsfs$count, S=jsfs$S, ts=ts, Ns=ns, cn=cn)
print(D)
print(jsfs)
#########################################################################################

# run model
mcmc.out <- ih_model_mcmc(par_init, params, mcmc.params, hp, D)
ll <- mcmc.out$mcmc.ll
mcmc.out <- mcmc.out$mcmc.out
#traj <- ih_model_plot2(mcmc.out, params)
#save(mcmc.out, traj, ll, file="mcmc_out_test.RData")


fout <- paste(outdir, "/intrahost.", inds, ".RData", sep="")
save(mcmc.out, ll, D, file=fout)

