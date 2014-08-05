args <- commandArgs(TRUE) ;
testno <- args[1]
expdir <- args[2]
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

expparamfile <- paste(expdir, "/exp.params", sep="")
paramfile <- paste(expdir, "/opt.params", sep="")
source(paramfile)
source(expparamfile)
allvars <- expand.grid(c(exp.params$dvars, exp.params$pvars))
varnames <- names(allvars)
	
switch(as.character(testno), 
	"1"={
		t_i <- seq(2,6,by=0.1)
		#for (v_it in seq(1, nrow(allvars)))
		#{
			source(paramfile)

		#	# output files
		#	varstr <- NULL
		#	for (i in seq(1, length(varnames)))
		#	{
		#		varstr <- paste(varstr, paste(varnames[i], "_", allvars[v_it, i], "_", sep=""), sep="")	
		#	}

		#	outdir <- paste(expdir, "/", varstr, sep="")
					

	
	
	
			llpdffile <- paste(outdir, "/ll_v_t.pdf", sep="")
			llfile <- paste(outdir, "/ll_v_t.dat", sep="")
			ll_all <- NULL
			for (j in exp.params$drep)
			{
				x <- read.dna(paste(outdir, "/seq.", j,".dat", sep=""))
				tstart <- 20
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
				
				#par_init <- ih_model_getparams(params.isopt, params.init)
				params <- c(params.init, list(nsites=nsites, Ns=ns))
				D <- list(mut=jsfs$mut, count=jsfs$count, S=jsfs$S, ts=ts, Ns=ns, cn=cn)
				print(D)
				print(jsfs)
				ll <- NULL
				for (t_ii in t_i)
				{
					params$t <- t_ii
	
					ll <- c(ll, ih_model_sfsll(params, D))
				}
				ll_all <- rbind(ll_all, ll)
			}
			write.table(ll_all, file=llfile, quote=FALSE, col.names=FALSE)
			q <- apply(ll_all, 2, quantile, c(0.25, 0.5, 0.75), na.rm=TRUE)
			pdf(llpdffile)
			plot(t_i, q[2,], type="l", lwd = 2, xlab="t", ylab="t_n")
			polygon(c(t_i, rev(t_i)), c(q[1,], rev(q[3,])), col="gray", border="gray")
			lines(t_i, q[2,], type="l", lwd = 2)
			dev.off()		
		#}
	}
)