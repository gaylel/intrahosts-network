#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

struct Itrajtype {
	double *i ; // fraction of population infected;
	double *s ;// fraction of population susceptible;
	double *t ; // times;
	double *inv_ne ; // inverse effective population size
	int I0 ;
	int N ; // length of I, t ;
	int length ;
	double T ; // length of infectious period ;
	double bet ;
	double gam ;
} ;


typedef struct Itrajtype Itrajtype;

struct ivectype {
	int * v ;
	int length ;
} ;

typedef struct ivectype ivectype ; 

/*
struct hosttype {
	phylo* tr ;
	int no ;
	int Nseq ;
	double tsamp ;
	double tinf ; // time of (first) infection
	double ts ;
	int *ti ;
	double *bt ;
} ;

typedef struct hosttype hosttype ;

hosttype* ih_hostinit(int ind, int Nseq, double tsamp) ;
hosttype* ih_hostinitts(hosttype* h, Itrajtype* I) ;
hosttype* ih_hostupdatets(hosttype* h, Itrajtype* I) ;
*/
Itrajtype * ih_drawI(double bet, double gam, int I0, int N, double ss ) ;
ivectype* ih_Ilookup(int l, Itrajtype* I) ;
int ih_tlookup(double t, Itrajtype* I) ;
double ih_interpI(double t, int t_i, int isI, Itrajtype* I) ;
double* ih_draw_ebt(Itrajtype* I, int n, int N, double t_from, double t_to) ;
double* ih_calculateess(double *ebt,  int n) ;
double ih_draw_nextbt(Itrajtype* I, int n, double t_from, double t_to) ;
