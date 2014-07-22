intrahost_R.so:	intrahost.c intrahost_R.c
	R CMD SHLIB intrahost_R.c intrahost.c

intrahost_R2.so:	intrahost2.c intrahost_R2.c sfs.c
	R CMD SHLIB intrahost_R2.c intrahost2.c sfs.c

rundir:
	mkdir ../run
	
#dstr := "date "+%Y%m%d"" \


runnew.%:	
	mkdir ../run/$*
	cat ../params/opt.params  > ../run/$*/opt.params
	
runall.%:
	cat ../submit.all.jobs |\
	 sed 's/rundir/$*/g' > ../run/$*/run.jobs

plotall.%:

