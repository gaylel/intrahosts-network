#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include "intrahost2.h"

//SEXP ih_R_drawSIR(SEXP R_bet, SEXP R_gam, SEXP R_I0, SEXP R_NS, SEXP R_ss, SEXP R_Nseq, SEXP R_ts, SEXP R_h) ;

SEXP ItrajtoR(Itrajtype *I) ;
sfsdatatype* ih_read_data(SEXP R_D) ;
SEXP ih_get_traj(SEXP R_bet, SEXP R_p, SEXP R_c, SEXP R_delt, SEXP R_V0, SEXP R_T0, SEXP R_ss) ;

sfsdatatype* ih_read_data(SEXP R_D)
{
	// assume R_D is in specific order ... can change this later
	SEXP R_mut, R_count, R_S, R_ts, R_ns, R_dim ;
	R_mut =  VECTOR_ELT(R_D, 0) ;
	R_count =  VECTOR_ELT(R_D, 1) ;
	R_S =  VECTOR_ELT(R_D, 2) ;
	R_ts =  VECTOR_ELT(R_D, 3) ;
	R_ns =  VECTOR_ELT(R_D, 4) ;
	R_dim = getAttrib(R_mut, R_DimSymbol) ;
	
	int** mut , *mut_c, lp, S, mut_n, *ns, i, j, *ptr;
	double *ts, *ptr2 ;
	sfsdatatype * x ;
	lp = INTEGER(coerceVector(R_dim, INTSXP))[0] ;
	mut_n = INTEGER(coerceVector(R_dim, INTSXP))[1] ;
	ptr = INTEGER(coerceVector(R_mut, INTSXP)) ;
	mut = calloc(mut_n, sizeof(int*)) ;
	for (i=0 ; i<mut_n ; i++)
	{
		mut[i] = calloc(lp, sizeof(int)) ;
		for (j=0 ; j<lp ; j++)
		{
			mut[i][j] = ptr[i*lp + j] ;
		}
	}
	
	ptr = INTEGER(coerceVector(R_count, INTSXP)) ;
	mut_c = calloc(mut_n, sizeof(int)) ;
	for (i=0 ; i<mut_n ; i++)
	{
		mut_c[i] = ptr[i] ;
	}
	
	ptr = INTEGER(coerceVector(R_ns, INTSXP)) ;
	ns = calloc(lp, sizeof(int)) ;
	for (i=0 ; i<lp ; i++)
	{
		ns[i] = ptr[i] ;
	}
	
	ptr2 = REAL(coerceVector(R_ts, REALSXP)) ;
	ts = calloc(lp, sizeof(double)) ;
	for (i=0 ; i<lp ; i++)
	{
		ts[i] = ptr2[i] ;
	}
	
	S = INTEGER(coerceVector(R_S, INTSXP))[0] ;
	
	x =  sfsdata_init(mut, mut_c, mut_n, lp, S, ts, ns) ;
	return x ;
}


SEXP ItrajtoR(Itrajtype *I)
{
	SEXP R_list, R_list_i, R_list_tc, R_list_v, R_list_ts, R_list_T, R_names ;
	
	PROTECT(R_list = allocVector(VECSXP, 5)) ;
	PROTECT(R_list_i = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_v = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_ts = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_tc = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_T = allocVector(REALSXP, 1)) ;
	
	
    

	double *i = REAL(R_list_i); // fraction of population infected;
	double *v = REAL(R_list_v);// fraction of population susceptible;
	double *tc = REAL(R_list_tc); // times;
	double *ts = REAL(R_list_ts); // times;
	double *T = REAL(R_list_T); // length of infectious period ;

	int j ;
	for (j=0 ; j<I->length ; j++)
	{
		i[j] = I->i[j] ;
		v[j] = I->v[j] ;
		tc[j] = I->tc[j] ;
		ts[j] = I->t[j] ;
	
	}
	T[0] = I->T ;
	SET_VECTOR_ELT(R_list, 0, R_list_v) ;
	SET_VECTOR_ELT(R_list, 1, R_list_tc) ;
	SET_VECTOR_ELT(R_list, 2, R_list_ts) ;
	SET_VECTOR_ELT(R_list, 3, R_list_i) ;
	SET_VECTOR_ELT(R_list, 4, R_list_T) ;
	
	PROTECT(R_names=allocVector(STRSXP,5)) ;
	SET_STRING_ELT(R_names, 0, mkChar("v")) ;
	SET_STRING_ELT(R_names, 1, mkChar("tc")) ;
	SET_STRING_ELT(R_names, 2, mkChar("t")) ;
	SET_STRING_ELT(R_names, 3, mkChar("i")) ;
	SET_STRING_ELT(R_names, 4, mkChar("T")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	UNPROTECT(7) ;
	return R_list ;
}

SEXP ih_get_traj(SEXP R_bet, SEXP R_p, SEXP R_c, SEXP R_delt, SEXP R_V0, SEXP R_T0, SEXP R_ss)
{
	double bet = REAL(coerceVector(R_bet, REALSXP))[0] ;
	double p = REAL(coerceVector(R_p, REALSXP))[0] ;
	double c = REAL(coerceVector(R_c, REALSXP))[0] ;
	double delt = REAL(coerceVector(R_delt, REALSXP))[0] ;
	double V0 = REAL(coerceVector(R_V0, REALSXP))[0] ;		// initial number of infected 
	double T0 = REAL(coerceVector(R_T0, REALSXP))[0] ;
	double ss = REAL(coerceVector(R_ss, REALSXP))[0] ;
	Itrajtype *I ;
	I = ih_drawI_Baccam(bet, delt, c, p, T0, V0, ss) ;
	SEXP R_I ;
	PROTECT(R_I = allocVector(VECSXP, 1)) ;
	SET_VECTOR_ELT(R_I, 0,  ItrajtoR(I)) ;
	ih_free(I) ;
	UNPROTECT(1) ; 
	return R_I ;
}


SEXP ih_R_getsfs(SEXP R_bet, SEXP R_p, SEXP R_c, SEXP R_delt, SEXP R_V0, SEXP R_T0, SEXP R_ss, SEXP R_t_off, SEXP R_D)
{
	/*******************************************************************************
	
	Draws deterministic SIR trajectory using Euler's forward method, and calculates 
	expected branching times, and site frequency spectrum 
	
	Args:	R_bet: growth rate
			R_gam: removal rate
			R-I0: initial number of infected
			R_NS:	Population size
			R_ss: stepsize
			R_Nseq: number of sequences sampled
			R_ts: time at which sequences are sampled, relative to infection at t=0
			
	Returns:
			pointer to Itrajtype 
	
	********************************************************************************/
	
	int i ;
	double bet = REAL(coerceVector(R_bet, REALSXP))[0] ;
	double p = REAL(coerceVector(R_p, REALSXP))[0] ;
	double c = REAL(coerceVector(R_c, REALSXP))[0] ;
	double delt = REAL(coerceVector(R_delt, REALSXP))[0] ;
	double V0 = REAL(coerceVector(R_V0, REALSXP))[0] ;		// initial number of infected 
	double T0 = REAL(coerceVector(R_T0, REALSXP))[0] ;
	double ss = REAL(coerceVector(R_ss, REALSXP))[0] ;
	sfsdatatype *sdt = ih_read_data(R_D) ;
	int Nsamp = sdt->lp ;
	double t_off = REAL(coerceVector(R_t_off, REALSXP))[0] ;
	Itrajtype* I ;
	int h_ind[] = {0, 1} ;
	int Niter = 5000 ;
	double *ll , S, *cn, *ptr;
	
	SEXP R_list, R_ll, R_names, R_pmean, R_flag, R_cn ;
	PROTECT(R_list = allocVector(VECSXP, 5) ) ;
	PROTECT(R_pmean = allocVector(REALSXP, 1)) ;
	PROTECT(R_flag = allocVector(INTSXP, 1)) ;
	PROTECT(R_ll = allocVector(REALSXP, 1)) ;
	PROTECT(R_cn = allocVector(REALSXP, Nsamp)) ;
	INTEGER(R_flag)[0] = 0 ;
	
	I = ih_drawI_Baccam(bet, delt, c, p, T0, V0, ss ) ;
	hosttype* h ;
	h = ih_hostinit(h_ind, sdt->ns, sdt->ts, Nsamp, t_off) ;
	
	// check params first
	cn = ih_get_cn(h, I) ;
	//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %i %8.4f %8.4f\n", h->t_off, h->t_inf, h->tsamp[0], h->tsamp[1], I->T, I->length, cn[0], cn[1]) ;
	if (ih_check(cn, h) > 0) 
	{
		ll = ih_psfs(sdt, I, h, Niter) ;
		//printf("ll = %8.4f\n, ", ll[0]) ;
		REAL(R_ll)[0] = ll[0] ;
		REAL(R_pmean)[0] = ll[1] ;
		ptr = REAL(R_cn) ;
		for (i=0 ; i<Nsamp ; i++)
		{
			ptr[i] = cn[i] ;
		}
		free(cn) ;
		free(ll) ;
		ih_hostfree(h) ;
		sfsdata_free(sdt) ;
		INTEGER(R_flag)[0] = 1 ;
	}
	SET_VECTOR_ELT(R_list, 0, R_ll) ;
	SET_VECTOR_ELT(R_list, 1,  ItrajtoR(I)) ;
	SET_VECTOR_ELT(R_list, 2,  R_pmean) ;
	SET_VECTOR_ELT(R_list, 3,  R_cn) ;
	SET_VECTOR_ELT(R_list, 4,  R_flag) ;
	PROTECT(R_names=allocVector(STRSXP,5)) ;
	SET_STRING_ELT(R_names, 0, mkChar("ll")) ;
	SET_STRING_ELT(R_names, 1, mkChar("sir")) ;
	SET_STRING_ELT(R_names, 2, mkChar("ES")) ;
	SET_STRING_ELT(R_names, 3, mkChar("cn")) ;
	SET_STRING_ELT(R_names, 4, mkChar("err")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	
	UNPROTECT(6) ;

	
	
	ih_free(I) ;
	return R_list ;
}





SEXP ih_R_drawtraj(SEXP R_bet, SEXP R_p, SEXP R_c, SEXP R_delt, SEXP R_V0, SEXP R_T0, SEXP R_ss, SEXP R_Nseq, SEXP R_ts, SEXP R_ti, SEXP R_h)
{
	/*******************************************************************************
	
	Draws deterministic SIR trajectory using Euler's forward method, and calculates 
	expected branching times, and site frequency spectrum 
	
	Args:	R_bet: growth rate
			R_gam: removal rate
			R-I0: initial number of infected
			R_NS:	Population size
			R_ss: stepsize
			R_Nseq: number of sequences sampled
			R_ts: time at which sequences are sampled, relative to infection at t=0
			
	Returns:
			pointer to Itrajtype 
	
	********************************************************************************/
	int i ;
	double bet = REAL(coerceVector(R_bet, REALSXP))[0] ;
	double p = REAL(coerceVector(R_p, REALSXP))[0] ;
	double c = REAL(coerceVector(R_c, REALSXP))[0] ;
	double delt = REAL(coerceVector(R_delt, REALSXP))[0] ;
	double V0 = REAL(coerceVector(R_V0, REALSXP))[0] ;		// initial number of infected 
	double T0 = REAL(coerceVector(R_T0, REALSXP))[0] ;
	double ss = REAL(coerceVector(R_ss, REALSXP))[0] ;
	int Nsamp = LENGTH(R_ts); 
	int *Nseq = calloc(Nsamp, sizeof(int)) ;
	double *ts = calloc(Nsamp, sizeof(double)) ; 
	double *ti = calloc(Nsamp, sizeof(double)) ; 
	
	for (i=0 ; i<Nsamp ; i++)
	{
		Nseq[i] = INTEGER(coerceVector(R_Nseq, INTSXP))[i] ;		
		ts[i] = REAL(coerceVector(R_ts, REALSXP))[i] ;
		ti[i] = REAL(coerceVector(R_ti, REALSXP))[i] ;
	}

	Itrajtype *I ; 
	ivectype *v ;
	int tsamp_i, Nav=1000 ;
	double *tsamp , *ess, pmean; 
	SEXP R_list, R_ess, R_ET, R_names, R_pmean, R_flag, R_cn ;
	PROTECT(R_list = allocVector(VECSXP, 6) ) ;
	PROTECT(R_ess = allocVector(REALSXP, Nseq[0] - 1)) ;
	PROTECT(R_ET = allocVector(REALSXP, Nseq[0] - 1)) ;
	PROTECT(R_pmean = allocVector(REALSXP, 1)) ;
	PROTECT(R_cn = allocVector(REALSXP, 1)) ;
	PROTECT(R_flag = allocVector(INTSXP, 1)) ;
	INTEGER(R_flag)[0] = 0 ;
	double *PR_ess = REAL(R_ess);
	double *PR_ET = REAL(R_ET);
	
	
	// draw trajectory
	I = ih_drawI_Baccam(bet, delt, c, p, T0, V0, ss ) ;
	hosttype* h ;
	h= ih_hostinit(0, Nseq, ts, Nsamp, ti[0]) ;
	
	
	
	
	

	
	
	
	int ind_from = ih_tlookup(ts[0], I) ;
	double ij = ih_interpI(ts[0], ind_from, 1, I) ;
	tsamp = ih_draw_ebt(I, Nseq[0], Nav, ts[0], 0) ;
	if (tsamp[0] >= 0)
	{
	
	ess = ih_calculateess(tsamp,  Nseq[0]) ;
	pmean = ih_nmutdist(tsamp, Nseq[0]) ;
	
	
	
	for (i=0 ; i<Nseq[0] - 1 ; i++)
	{
		PR_ess[i] = ess[i] ;
		PR_ET[i] = tsamp[i] ;
	}
	REAL(R_pmean)[0] = pmean ;
	REAL(R_cn)[0] = ij ;
	INTEGER(R_flag)[0] = 1 ;
	free(ess) ;
	free(tsamp) ;
	}
	
	SET_VECTOR_ELT(R_list, 0, R_ess) ;
	SET_VECTOR_ELT(R_list, 1,  ItrajtoR(I)) ;
	SET_VECTOR_ELT(R_list, 2,  R_ET) ;
	SET_VECTOR_ELT(R_list, 3,  R_pmean) ;
	SET_VECTOR_ELT(R_list, 4,  R_cn) ;
	SET_VECTOR_ELT(R_list, 5,  R_flag) ;
	PROTECT(R_names=allocVector(STRSXP,6)) ;
	SET_STRING_ELT(R_names, 0, mkChar("ESS")) ;
	SET_STRING_ELT(R_names, 1, mkChar("sir")) ;
	SET_STRING_ELT(R_names, 2, mkChar("ET")) ;
	SET_STRING_ELT(R_names, 3, mkChar("ES")) ;
	SET_STRING_ELT(R_names, 4, mkChar("cn")) ;
	SET_STRING_ELT(R_names, 5, mkChar("err")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	
	UNPROTECT(7) ;
	
	ih_free(I) ;
	free(Nseq) ;
	free(ts) ;

	return R_list ; 
	
}



SEXP ih_R_drawtraj2(SEXP R_bet, SEXP R_p, SEXP R_c, SEXP R_delt, SEXP R_V0, SEXP R_T0, SEXP R_ss, SEXP R_Nseq, SEXP R_ts, SEXP R_ti, SEXP R_h)
{
	/*******************************************************************************
	
	Draws deterministic SIR trajectory using Euler's forward method, and calculates 
	expected branching times, and site frequency spectrum 
	
	Args:	R_bet: growth rate
			R_gam: removal rate
			R-I0: initial number of infected
			R_NS:	Population size
			R_ss: stepsize
			R_Nseq: number of sequences sampled
			R_ti: time at which sequences are sampled, relative to infection at t=0
			R_ts: absolute times of sampling
			
	Returns:
			pointer to Itrajtype 
	
	********************************************************************************/
	int i ;
	srand(time(NULL)) ;
	double bet = REAL(coerceVector(R_bet, REALSXP))[0] ;
	double p = REAL(coerceVector(R_p, REALSXP))[0] ;
	double c = REAL(coerceVector(R_c, REALSXP))[0] ;
	double delt = REAL(coerceVector(R_delt, REALSXP))[0] ;
	double V0 = REAL(coerceVector(R_V0, REALSXP))[0] ;		// initial number of infected 
	double T0 = REAL(coerceVector(R_T0, REALSXP))[0] ;
	double ss = REAL(coerceVector(R_ss, REALSXP))[0] ;
	int Nsamp = LENGTH(R_ts); 
	int *Nseq = calloc(Nsamp, sizeof(int)) ;
	double *ts = calloc(Nsamp, sizeof(double)) ; 
	double ti = REAL(coerceVector(R_ti, REALSXP))[0] ;
	
	for (i=0 ; i<Nsamp ; i++)
	{
		Nseq[i] = INTEGER(coerceVector(R_Nseq, INTSXP))[i] ;		
		ts[i] = REAL(coerceVector(R_ts, REALSXP))[i] ;
		
	}

	Itrajtype *I ; 
	ivectype *v ;
	int tsamp_i, Nav=1000 ;
	double *tsamp , *ess, pmean; 
	SEXP R_list, R_ess, R_ET, R_names, R_pmean, R_flag, R_cn ;
	PROTECT(R_list = allocVector(VECSXP, 6) ) ;
	PROTECT(R_ess = allocVector(REALSXP, Nseq[0] - 1)) ;
	PROTECT(R_ET = allocVector(REALSXP, Nseq[0] - 1)) ;
	PROTECT(R_pmean = allocVector(REALSXP, 1)) ;
	PROTECT(R_cn = allocVector(REALSXP, 1)) ;
	PROTECT(R_flag = allocVector(INTSXP, 1)) ;
	INTEGER(R_flag)[0] = 0 ;
	double *PR_ess = REAL(R_ess);
	double *PR_ET = REAL(R_ET);
	
	
	// draw trajectory
	I = ih_drawI_Baccam(bet, delt, c, p, T0, V0, ss ) ;
	hosttype* h ;
	int h_ind[] = {0, 1} ;
	h= ih_hostinit(h_ind, Nseq, ts, Nsamp, ti) ;
	//h = ih_draw_coaltree(I, h, 100) ;
	
	
	
	

	
	
	
	int ind_from = ih_tlookup(ts[0], I) ;
	double ij = ih_interpI(ts[0], ind_from, 1, I) ;
	tsamp = ih_draw_ebt(I, Nseq[0], Nav, ts[0], 0) ;
	if (tsamp[0] >= 0)
	{
	
	ess = ih_calculateess(tsamp,  Nseq[0]) ;
	pmean = ih_nmutdist(tsamp, Nseq[0]) ;
	
	
	
	for (i=0 ; i<Nseq[0] - 1 ; i++)
	{
		PR_ess[i] = ess[i] ;
		PR_ET[i] = tsamp[i] ;
	}
	REAL(R_pmean)[0] = pmean ;
	REAL(R_cn)[0] = ij ;
	INTEGER(R_flag)[0] = 1 ;
	free(ess) ;
	free(tsamp) ;
	}
	
	SET_VECTOR_ELT(R_list, 0, R_ess) ;
	SET_VECTOR_ELT(R_list, 1,  ItrajtoR(I)) ;
	SET_VECTOR_ELT(R_list, 2,  R_ET) ;
	SET_VECTOR_ELT(R_list, 3,  R_pmean) ;
	SET_VECTOR_ELT(R_list, 4,  R_cn) ;
	SET_VECTOR_ELT(R_list, 5,  R_flag) ;
	PROTECT(R_names=allocVector(STRSXP,6)) ;
	SET_STRING_ELT(R_names, 0, mkChar("ESS")) ;
	SET_STRING_ELT(R_names, 1, mkChar("sir")) ;
	SET_STRING_ELT(R_names, 2, mkChar("ET")) ;
	SET_STRING_ELT(R_names, 3, mkChar("ES")) ;
	SET_STRING_ELT(R_names, 4, mkChar("cn")) ;
	SET_STRING_ELT(R_names, 5, mkChar("err")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	
	UNPROTECT(7) ;
	
	ih_free(I) ;
	free(Nseq) ;
	free(ts) ;

	return R_list ; 
	
}



