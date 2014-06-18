#	i/o functions 
#
#

ih_io_readdata <- function()
{
	# read in data
	data(HorseFlu)
	x <- HorseFlu
	return(x)
}


ih_io_sfs <- function(x, uniqseq, cseq, inds)
{
		cs <- as.character(uniqseq@uniqdna[match(cseq, labels(uniqseq@uniqdna)),])
		xm <- as.matrix(x@dna@dna[[1]])
		bf <- numeric(nsites)
		ids <- ih_io_getseqIDs(x, inds)
		nseq <- length(ids)
		nseq_all <- nrow(xm)
		print(ids)
		for (i in seq(1, nsites))
		{
			bf[i] <- as.integer((1 - base.freq(xm[ids,i])[[cs[,i]]])*nseq)
		}
		print(bf)
		sfs <- tabulate(bf, nbins= nseq_all - 1)
		#sfs <- tabulate(bf)
		
		return(sfs)
}

ih_io_getseqIDs <- function(x, inds)
{
	ids <- subset(x, individual=inds)
	return(row.names(ids@dna@meta))
}
