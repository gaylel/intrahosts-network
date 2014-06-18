library("OutbreakTools")
library("ape")

#input/output functions 
dyn.load("intrahost_R.so") 
source("ih_io.R")
source("ih_model.R")


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
inds <- c(16)
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

t <- 20
params <- list(bet= 2, gam= 0.5, I0=10, NS=100000)
ih <- ih_model_fit(t, ns[1], params)
ll <- ih_model_ll(ih, sfs[1:(ns[1]-1)])
print(ll)
ih_model_plot(ih, t)

