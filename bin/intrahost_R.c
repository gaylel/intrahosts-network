#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include "intrahost.h"

SEXP ih_R_drawSIR(SEXP R_bet, SEXP R_gam, SEXP R_I0, SEXP R_NS, SEXP R_ss, SEXP R_Nseq, SEXP R_ts) ;
SEXP ItrajtoR(Itrajtype *I) ;

SEXP ItrajtoR(Itrajtype *I)
{
	SEXP R_list, R_list_i, R_list_s, R_list_t, R_list_T, R_names ;
	
	PROTECT(R_list = allocVector(VECSXP, 4)) ;
	PROTECT(R_list_i = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_s = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_t = allocVector(REALSXP, I->length)) ;
	PROTECT(R_list_T = allocVector(REALSXP, 1)) ;
	
	
    

	double *i = REAL(R_list_i); // fraction of population infected;
	double *s = REAL(R_list_s);// fraction of population susceptible;
	double *t = REAL(R_list_t); // times;
	double *T = REAL(R_list_T); // length of infectious period ;

	int j ;
	for (j=0 ; j<I->length ; j++)
	{
		i[j] = I->i[j] ;
		s[j] = I->s[j] ;
		t[j] = I->t[j] ;
	}
	T[0] = I->T ;
	SET_VECTOR_ELT(R_list, 0, R_list_i) ;
	SET_VECTOR_ELT(R_list, 1, R_list_s) ;
	SET_VECTOR_ELT(R_list, 2, R_list_t) ;
	SET_VECTOR_ELT(R_list, 3, R_list_T) ;
	
	PROTECT(R_names=allocVector(STRSXP,4)) ;
	SET_STRING_ELT(R_names, 0, mkChar("i")) ;
	SET_STRING_ELT(R_names, 1, mkChar("s")) ;
	SET_STRING_ELT(R_names, 2, mkChar("t")) ;
	SET_STRING_ELT(R_names, 3, mkChar("T")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	UNPROTECT(6) ;
	return R_list ;
}


SEXP ih_R_drawSIR(SEXP R_bet, SEXP R_gam, SEXP R_I0, SEXP R_NS, SEXP R_ss, SEXP R_Nseq, SEXP R_ts)
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

	double bet = REAL(coerceVector(R_bet, REALSXP))[0] ;
	double gam = REAL(coerceVector(R_gam, REALSXP))[0] ;
	int I0 = INTEGER(coerceVector(R_I0, INTSXP))[0] ;		// initial number of infected 
	int NS = INTEGER(coerceVector(R_NS, INTSXP))[0] ;		// size of intrahost population
	double ss = REAL(coerceVector(R_ss, REALSXP))[0] ;
	int Nseq = INTEGER(coerceVector(R_Nseq, INTSXP))[0] ;		
	double ts = REAL(coerceVector(R_ts, REALSXP))[0] ;
	

	Itrajtype *I ; 
	ivectype *v ;
	int tsamp_i, Nav=1000 ;
	double *tsamp , *ess; 
	
	// draw SIR curve
	I = ih_drawI(bet, gam, I0, NS, ss) ; 
	tsamp = ih_draw_ebt(I, Nseq, Nav, ts, 0) ;
	ess = ih_calculateess(tsamp,  Nseq) ;
	
	SEXP R_list, R_ess, R_ET, R_names ;
	PROTECT(R_list = allocVector(VECSXP, 3) ) ;
	PROTECT(R_ess = allocVector(REALSXP, Nseq - 1)) ;
	PROTECT(R_ET = allocVector(REALSXP, Nseq - 1)) ;
	double *PR_ess = REAL(R_ess);
	double *PR_ET = REAL(R_ET);
	int i ;
	for (i=0 ; i<Nseq - 1 ; i++)
	{
		PR_ess[i] = ess[i] ;
		PR_ET[i] = tsamp[i] ;
	}
	SET_VECTOR_ELT(R_list, 0, R_ess) ;
	SET_VECTOR_ELT(R_list, 1,  ItrajtoR(I)) ;
	SET_VECTOR_ELT(R_list, 2,  R_ET) ;
	
	PROTECT(R_names=allocVector(STRSXP,3)) ;
	SET_STRING_ELT(R_names, 0, mkChar("ESS")) ;
	SET_STRING_ELT(R_names, 1, mkChar("sir")) ;
	SET_STRING_ELT(R_names, 2, mkChar("ET")) ;
	setAttrib(R_list, R_NamesSymbol, R_names);
	
	UNPROTECT(4) ;
	return R_list ; 
	
}

