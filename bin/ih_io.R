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

ih_io_readcndata <- function(fname)
{
	cn <- read.csv(fname)
	return(cn)
}

ih_io_convertdates <- function(sdates, tstart)
{
	rv <- as.numeric(sdates - sdates[1])
	return(rv)
}


ih_io_sfs <- function(x, uniqseq, cseq, inds)
{
		cs <- as.character(uniqseq@uniqdna[match(cseq, labels(uniqseq@uniqdna)),])
		xm <- as.matrix(x@dna@dna[[1]])
		nsites <- length(get.dna(x)[[1]][1,])

		bf <- numeric(nsites)
		ids <- ih_io_getseqIDs(x, inds)
		nseq <- length(ids)
		nseq_all <- nrow(xm)
		print(ids)
		for (i in seq(1, nsites))
		{
			bf[i] <- as.integer(round((1 - base.freq(xm[ids,i])[[cs[,i]]])*nseq))
		}
		print(bf)
		sfs <- tabulate(bf, nbins= nseq_all - 1)
		#sfs <- tabulate(bf)
		
		return(sfs)
}

ih_io_jointsfs <- function(x, uniqseq, cseq, inds)
{
	cs <- as.character(uniqseq@uniqdna[match(cseq, labels(uniqseq@uniqdna)),])
	xm <- as.matrix(x@dna@dna[[1]])
	bf <- matrix(0, length(inds), nsites)
	nsites <- length(get.dna(x)[[1]][1,])
	for (j in seq(1,length(inds)))
	{
		ids <- ih_io_getseqIDs(x, inds[j])
		nseq <- length(ids)
		for (i in seq(1, nsites))
		{
			bf[j, i] <- as.integer(round((1 - base.freq(xm[ids,i])[[cs[,i]]])*nseq))
		}
	}
	return(bf)
}

ih_get_consensus<- function(x, uniqseq, sids)
{
	# get ids for first time sample
	ids <- row.names(x@dna@meta[which(x@dna@meta[,"sampleID"] == sids[1]),])
	uu1 <- match(ids, unlist(uniqseq@uniqID))
	uid <- unlist(uniqseq@uniqID)
	uidc <- do.call(rbind, lapply(uniqseq@uniqID, function(x) length(x)))
	names(uid) <- rep(names(uniqseq@uniqID),uidc)
	us1 <- uid[uu1]
	
	# frequencies of unique sequences
	f1 <- rle(sort(names(us1)))
	uus1 <- uniqseq@uniqdna[unique(names(us1)),]
	
	lab11 <- f1$lengths[match(unique(names(us1)), f1$values)]
	m <- which(lab11==max(lab11))		
	cseq <-  row.names(uus1)[m]
	return(cseq)		
}

ih_io_jsfs <- function(x, uniqseq, cseq, sids)
{
	cs <- as.character(uniqseq@uniqdna[match(cseq, labels(uniqseq@uniqdna)),])
	xm <- as.matrix(x@dna@dna[[1]])
	bf <- matrix(0, length(sids), nsites)
	nsites <- length(get.dna(x)[[1]][1,])
	for (j in seq(1,length(sids)))
	{
		ids <- row.names(x@dna@meta[which(x@dna@meta[,"sampleID"] == sids[j]),])
		nseq <- length(ids)
		for (i in seq(1, nsites))
		{
			bf[j, i] <- as.integer(round((1 - base.freq(xm[ids,i])[[cs[,i]]])*nseq))
		}
	}
	
	ujsfs <- unique(bf, MARGIN=2)
	ujsfs <- matrix(ujsfs[, -which(colSums(ujsfs)==0)],nrow=length(sids))
	ujsfs_c <- tabulate(match(data.frame(bf), data.frame(ujsfs)))
	S = sum(colSums(ujsfs) * ujsfs_c)
	jsfs <- list(mut=ujsfs, count=ujsfs_c, S=S)
	return(jsfs)
}


ih_network <- function(bf)
{
	nw <- matrix(0, nrow(bf), nrow(bf))
	for (j in seq(1, ncol(bf)))
	{
		ids <- which(bf[, j]>=1)
		if (length(ids) > 1)
		{
			for (m in seq(1, length(ids)-1))
			{
				for (n in seq(m+1, length(ids)))
				{
					nw[ids[m], ids[n]] <- 1
				}
			}
		}
	}
	return(nw)
}

ih_network2 <- function(bf)
{
	# plot 
}



ih_io_getseqIDs <- function(x, inds)
{
	ids <- subset(x, individual=inds)
	return(row.names(ids@dna@meta))
}
