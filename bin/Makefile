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


testdir:
	mkdir ../test

testnew.%:	
	mkdir ../test/$*
	cat ../params/data.params > ../test/$*/data.params
	cat ../params/opt.params  > ../test/$*/opt.params

testdata.%:	
	Rscript intrahostdata.R ../test/$*/data.params ../test/$*
	#seq-gen -mHKY -s7e-5 -l12000 < ../test/$*/tr.nwk > ../test/$*/seq.dat
	./seqgendata.sh ../test/$*/tr.nwk ../test/$*/seq.dat ../test/$*/data.params 
	#echo " 100 903" > ../test/$*/seq.dat
	#cat ../test/$*/seq.dat.all | grep ^t >> ../test/$*/seq.dat
	#cat ../test/$*/seq.dat.all | grep -v ^t > ../test/$*/nodes.dat
	
	
testrun.%:	 
	Rscript intrahosttest.R ../test/$*/opt.params ../test/$*/ ../test/$*/

testall.%:
	cat ../submit.test.jobs |\
	sed 's/rundir/$*/g' > ../test/$*/run.jobs

expnew.%:
	mkdir ../test/$*
	cat ../params/data.params > ../test/$*/data.params
	cat ../params/opt.params  > ../test/$*/opt.params
	cat ../params/exp.params > ../test/$*/exp.params
	


