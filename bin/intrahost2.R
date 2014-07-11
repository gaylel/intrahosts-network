library("OutbreakTools")
library("ape")
library("MCMCpack")

#input/output functions 
dyn.load("intrahost_R2.so") 
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
cseq <- row.names(IDcounts)[2]

#inds <- get.individuals(x)


inds <- c(25)
#inds <- seq(1,50)
sfs <- ih_io_sfs(x, uniqseq, cseq, inds)
#print(sfs)

seqs.meta <- subset(x, individual=inds)@dna@meta	
seqs.sids <- unique(seqs.meta$sampleID)
seqs.dates <- seqs.meta[match(seqs.sids, seqs.meta$sampleID), "date"]
o <- order(seqs.dates)
seqs.dates <- seqs.dates[o]
seqs.sids <- seqs.sids[o]
ns <- numeric(length(seqs.sids)) 
for (i in seq(1, length(seqs.sids)))
{
	ns[i] <- length(which(seqs.meta[,"sampleID"]==seqs.sids[i]))
}

jsfs <- ih_io_jsfs(x, uniqseq, cseq, seqs.sids)



print(seqs.dates)
print(jsfs)

t <- 3
S <- sum(sfs)

print(S)


bet = 3.4e-5
p=7.9e-3
c=3.3
delt=3.4
V0 = 3.5e-1
T0 = 4e8
mr=1e-4

print(ns)
par_init <-c(t,bet,c,p,delt,mr) 
params <- list(bet= bet, p=p, c=c, delt=delt, V0=V0, T0=T0,  nsites=nsites, Ns = ns, t=t, mr=mr)
V <- diag(length(par_init))
mcmc.params <- list(tune=c(0.05, 0.02, 0.02, 0.02, 0.02, 0.01), V=V, verbose=100, burnin=10000, mcmc=50000)
#data <- list(sfs=sfs[1:(ns[1]-1)], S=S)
D <- list(mut=jsfs$mut, count=jsfs$count, S=jsfs$S, ts=c(20, 22), Ns=ns)
#is.list(D)

#print(D)
#print(str(D))
#ih_model_sfsll(params, D)

mcmc.out <- ih_model_mcmc(par_init, params, mcmc.params, D)
traj <- ih_model_plot2(mcmc.out, params)
save(mcmc.out, traj, file="mcmc_out12.RData")

