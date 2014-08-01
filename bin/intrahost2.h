#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "sfs.h"

struct Itrajtype {
	double *tc ; // number of uninfected target cells ;
	double *v ; // infectious viral titer ;
	double *i ; // number of productively infected cells ;
	double *t ; // times
	double p ;	// shedding rate
	double c ;	// clearance rate
	double bet ;	//	cell infection rate
	double delt ; 	// 	cell death rate
	double *inv_ne ; // inverse effective population size
	int T0 ; 
	int V0 ;
	int length ;
	int N ; // length of I, t ;
	double T ; // length of infectious period ;
} ;

typedef struct Itrajtype Itrajtype ;


struct ivectype {
	int * v ;
	int length ;
} ;

typedef struct ivectype ivectype ; 


struct hosttype {
	int *no ;	// host number
	int Nsamp ; 	// number of samples taken
	int *Nseq ;	//	number of sequences for each sample
	double *tsamp ;	// sampling times
	double t_inf ; // time of (first) infection
	double t_off ;	// time to sampling
} ;

typedef struct hosttype hosttype ;

struct bttype {
	int *nlin ;	// number of lineages
	int *h1 ; // host index
	int *h2 ;
	double *t ;	// branching time 
	int *k ; // key 1=coalescent, -n=adding n lineages
	int N ;
} ;

typedef struct bttype bttype ; 

struct sfsdatatype {
	int **mut ;		// mutation patterns
	int *mutc ; 	// counts
	int mut_n ; 	// number of patterns
	int lp ; 		// length of pattern
	int S ; 		// number of mutations
	double *ts ; 	// absolute times of sampling
	int* ns	;		// number of sequences at each sampling point 
} ;

typedef struct sfsdatatype sfsdatatype ;

hosttype* ih_hostinit(int* ind, int* Nseq, double *tsamp, int Nsamp, double t_off) ;
/*
hosttype* ih_hostinitts(hosttype* h, Itrajtype* I) ;
hosttype* ih_hostupdatets(hosttype* h, Itrajtype* I) ;
*/
ivectype* ih_Ilookup(int l, Itrajtype* I) ;
int ih_tlookup(double t, Itrajtype* I) ;
double ih_interpI(double t, int t_i, int isI, Itrajtype* I) ;
double* ih_draw_ebt(Itrajtype* I, int n, int N, double t_from, double t_to) ;
double* ih_calculateess(double *ebt,  int n) ;
double ih_draw_nextbt(Itrajtype* I, int n, double t_from, double t_to) ;
double ih_nmutdist(double *ebt, int n) ;
void ih_free(Itrajtype* I) ;
void ih_hostfree(hosttype *h) ;
Itrajtype * ih_drawI_Baccam(double bet, double delt, double c, double p, double T0, double V0, double ss ) ;
Itrajtype * ih_drawI_Baccam2(double bet, double delt, double c, double p, double T0, double V0, double ss ) ;
sfstype* ih_draw_coaltree(Itrajtype *I, hosttype* h, int Niter) ;
sfstype* ih_draw_coaltree2(Itrajtype *I, hosttype* h, int Niter, sfsdatatype *sdt) ;
sfsdatatype* sfsdata_init(int** mut, int *mutc, int mut_n, int lp, int S, double* ts, int *ns) ;
bttype * ih_btmxreset(bttype *x) ;
bttype * ih_btmxupdate2(bttype* bt, int nlin, int h1, int h2, double t, int k) ;
void ih_btfree(bttype *x) ;
double* sfs_calcll(sfstype* x, sfsdatatype *d, int Niter) ;
double* sfs_calcll2(sfstype* x, sfsdatatype *d, int Niter) ;
void sfsdata_free(sfsdatatype *x ) ;
double* ih_psfs(sfsdatatype* sdt, Itrajtype *I, hosttype *h, int Niter) ;
double* ih_psfs2(sfsdatatype* d, Itrajtype *I, hosttype *h, int Niter) ;
double* ih_get_cn(hosttype* h, Itrajtype *I) ;
double* ih_get_ebt(sfsdatatype* sdt, Itrajtype *I, hosttype *h, int Niter) ;
sfstype* ih_draw_allbt(int Niter, hosttype* h, Itrajtype* I) ;
int ih_check(double* cn, hosttype *h) ;
double ih_calculateess1(double *ebt, int n) ;
double ih_calculateess2(double *ebt, int n, int i, double sumkET) ;