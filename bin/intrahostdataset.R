## generate some data

args <- commandArgs(TRUE) ;
expdir <- args[1]

library("OutbreakTools")
library("ape")

#input/output functions 
dyn.load("intrahost_R2.so") 
source("ih_io.R")
source("ih_model3.R")

paramfile <- paste(expdir, "/data.params", sep="")


expparamfile <- paste(expdir, "/exp.params", sep="")
source(expparamfile)

set.seed(exp.params$seed)
allvars <- expand.grid(c(exp.params$dvars, exp.params$pvars))
varnames <- names(allvars)
	
v_it <- 1
for (v_it in seq(1, nrow(allvars)))
{
	source(paramfile)

	# output files
	varstr <- NULL
	for (i in seq(1, length(varnames)))
	{
		varstr <- paste(varstr, paste(varnames[i], "_", allvars[v_it, i], "_", sep=""), sep="")	
	}

	outdir <- paste(expdir, "/", varstr, sep="")
	dir.create(outdir)

	# create /data.params file
	dparamfile <- paste(outdir, "/data.params", sep="") 
	sysstr <- paste("cat", paramfile)
	for (i in seq(1, length(varnames)))
	{
		sysstr <- paste(sysstr, paste("| sed 's/", varnames[i], ".*/", varnames[i], "=", allvars[v_it, i], ",/'", sep=""), sep="")
	}
	sysstr <- paste(sysstr, paste( "> ", dparamfile), sep="")
	system(sysstr) 

	tstart <- 20
	for (i in seq(1, length(varnames)))
	{
		if (varnames[i] %in% exp.params$pvars)
		{
			params[varnames[i]] <- allvars[v_it, i]
		}
		
		if (varnames[i] %in% exp.params$dvars)
		{
			D[varnames[i]] <- allvars[v_it, i]
		}
	}


	# sample trajectory
	params$bet <- params$bet/params$sc
	params$p <- params$p * params$sc
	
	for (j in exp.params$drep)
	{
		
		trfile <- paste(outdir, "/tr.", j,".nwk", sep="")
		datafile <- paste(outdir, "/dataparams.", j, ".RData", sep="")
		seqfile <- paste(outdir, "/seq.", j, ".dat", sep="")

		
		sir <- ih_get_traj_R(params$bet, params$p, params$c, params$delt, params$V0, params$T0, params$ss)
		sir <- sir[[1]]

		# sample branching times

		ebt <- ih_model_get_ebt(params, D, 1)
		print(ebt)
		# create phylo from branching times
		tr <- ih_model_get_phylo((rev(ebt)), D$Ns, D$ts)

		print(tr)
		write.tree(tr, file=trfile)

		save(sir, ebt, tr, params, D, file=datafile)


		sysstr <- paste("./seqgendata.sh ", trfile, seqfile, dparamfile)
		print(sysstr)
		system(sysstr)
	}
}