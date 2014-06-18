#include "intrahost.h"
#include <Rmath.h>
/*
hosttype* ih_hostinit(int ind, int Nseq, double tsamp)
{
	hosttype *h ;
	h = calloc(1, sizeof(hosttype)) ;
	
	h->no = ind ;
	h->Nseq = Nseq ;
	h->tsamp = tsamp ;
	return h ;
}

hosttype* ih_hostinitts(hosttype* h, Itrajtype* I)
{
	ivectype *v ;
	//int tsamp_i ;
	v = ih_Ilookup(h->Nseq, I) ; // indices get time interval (relative to trajectory start) in which sample is taken.
	//tsamp_i = (v->length > 0) ? v->v[v->length - 1] : -1 ;
	h->ti = calloc(2, sizeof(int)) ;
	h->ti[0] = v->v[0] ;
	h->ti[1] = v->v[1] ;  
	return h ;
}

hosttype* ih_hostupdatets(hosttype* h, Itrajtype* I) 
{
	//h->ts = I->t[tsamp_i] ;
	h->ts =20 ;
	return h ;
}
*/

Itrajtype * ih_drawI(double bet, double gam, int I0, int N, double ss )
{
	/*******************************************************************************
	
	Draws deterministic SIR trajectory using Euler's forward method
	
	Args:	beta: growth rate
			gamma: removal rate
			I0: initial number of infected
			N:	Population size
			ss: stepsize
			
	Returns:
			pointer to Itrajtype 
	
	********************************************************************************/
	
	double thr = (double ) 1 / N ;
	double t_i = 0, d_s, d_i;
	int n=0;
	
	Itrajtype * I ;
	I = calloc(1, sizeof(Itrajtype));
	I->length = 1 ;
	I->I0 = 1 ;
	I->i = calloc(1, sizeof(double)) ; 
	I->i[0] = (double ) I0 / N ;
	I->s = calloc(1, sizeof(double)) ; 
	I->s[0] = 1.0 - I->i[0];
	I->t = calloc(1, sizeof(double)) ;
	I->t[0] = 0 ;
	I->inv_ne = calloc(1, sizeof(double)) ; 
	I->inv_ne[0] = I->i[0] / (2 * bet * I->s[0]) ;
	I->N = N ;
	I->bet = bet ; 
	I->gam = gam ;
	
	
	while (I->i[n] >= thr)
	{
		//printf("%8.4f %8.4f %8.4f\n", I->t[n], I->i[n] * I->N, 1/I->inv_ne[n]) ;
		d_s = -(bet * I->s[n] * I->i[n]) * ss  ;
		d_i = -d_s - ((gam * I->i[n]) * ss) ;
		I->i = realloc(I->i,  (I->length + 1) * sizeof(double)) ;
		I->s = realloc(I->s,  (I->length + 1) * sizeof(double)) ;
		I->t = realloc(I->t,  (I->length + 1) * sizeof(double)) ;
		I->inv_ne = realloc(I->inv_ne,  (I->length + 1) * sizeof(double)) ;
		I->length++ ;
		n++ ;
		I->i[n] = I->i[n-1] + d_i ;
		I->t[n] = I->t[n-1] + ss ;
		I->s[n] = I->s[n-1] + d_s ; 
		I->inv_ne[n] = (2 * bet * I->s[n]) / ((I->i[n] * I->N) - 1);
		I->T+= ss ; 
		 
	}
	
	return I ;
		
}

ivectype* ih_Ilookup(int l, Itrajtype* I)
{
/*******************************************************************************
	
	Gets index of deterministic SIR trajectory for value of I*N
	
	Args:	l number of sequences
			I Itraj
			
	Returns:
			index of I 
	
	********************************************************************************/
	int n, nfound=0, *ind ;
	ivectype* v = calloc(1, sizeof(ivectype)) ;
	for (n=1 ; n< (I->length - 1) ; n++)
	{
		if ((((I->N * I->i[n]) > (double) l) && ((I->N * I->i[n-1]) <= (double) l)) || (((I->N * I->i[n+1]) <= (double) l) && ((I->N * I->i[n]) > (double) l)))
		{
			//printf("%8.4f, %8.4f\n", I->t[n], I->i[n]) ;
			if (nfound == 0)
			{
				v->v = calloc(1, sizeof(int)) ;
			}
			else
			{
				v->v = realloc(v->v, (nfound+1)*sizeof(int)) ;
				
			}
			v->v[nfound] = n ;
			//printf("ind %i, %8.4f, %8.4f\n", n, I->t[n], I->i[n]) ;
			nfound++ ;
			v->length = nfound ;
		}
	}
	return v ;
}

int ih_tlookup(double t, Itrajtype* I)
{
/*******************************************************************************
	
	Gets index of deterministic SIR trajectory for value of t
	
	Args:	t time to look up
			I Itraj
			
	Returns:
			index of I (lower bound) 
	
	********************************************************************************/
	double ss = I->t[1] - I->t[0] ;
	int ind = (int) floor(t / ss) ;
	return ind ;
}

double ih_interpI(double t, int t_i, int isI, Itrajtype* I)
{
	double diff, frac_t, rv ;
	frac_t = I->t[t_i+1] - I->t[t_i] ;
	frac_t = (t - I->t[t_i]) / frac_t ;   
	switch(isI)
	{
		case 0: // S
			diff = I->s[t_i + 1] - I->s[t_i] ;
			rv = I->s[t_i] + frac_t * diff ; 	 
		break ;
		case 1: // I
			diff = I->i[t_i + 1] - I->i[t_i] ;
			rv = I->i[t_i] + frac_t * diff ;
		break ;
	}
	return rv ;
}

double* ih_draw_ebt(Itrajtype* I, int n, int N, double t_from, double t_to)
{
	double *bt =calloc(N, sizeof(double)) , *ebt = calloc(n-1, sizeof(double));
	int ni = n, i, j;
	double last = t_from ; 
	printf("\n\n") ; 
	for (j=n-2 ; j>=0 ; j--)
	{
		ebt[j] = 0 ;
	for (i=0 ; i<N ; i++)
	{
		bt[i] = ih_draw_nextbt(I, j + 2, last, t_to) ;
		ebt[j] += bt[i] ;
		//printf("%8.4f\n", bt[i]) ;
	}
		ebt[j] = ebt[j] / N ;
		ebt[j] = last - ebt[j] ;
		last = last - ebt[j] ; 
		 //printf("%8.4f\t%8.4f\n", ebt[j], last) ;
	}
	return ebt;
}

double* ih_calculateess(double *ebt,  int n)
{
	// calculate expected site spectrum
	int L = n -1, i, k;
	double sumkET = 0 ;
	double* ess = calloc(L, sizeof(double)) ;
	for (i=0 ; i<L ; i++)
	{
		sumkET += (i+2) * ebt[i] ; 
	}	
	
	for (i=1 ; i<=L ; i++)
	{
		ess[i-1] = 0 ;
		for (k=2 ; k<= (n - i + 1) ; k++)
		{
			ess[i-1] += k * (k-1) *choose(n-k, i-1) * ebt[k-2] ;
		}
		ess[i-1] = ess[i-1] * beta(n-i,i) / sumkET ;
		//printf("%8.4f, ", ess[i-1]) ; 
	}
	return ess ;
}

double ih_draw_nextbt(Itrajtype* I, int n, double t_from, double t_to)
{
	// draw expected coalescence time
	int ind_from, ind_to ;
	int nchoose2 = (n * (n-1)) / 2 ;
	ind_from = ih_tlookup(t_from, I) ;
	double sj = ih_interpI(t_from, ind_from, 0, I) ;
	double ij = ih_interpI(t_from, ind_from, 1, I) ;
	ind_to = ih_tlookup(t_to, I) ;
	double si = ih_interpI(t_to, ind_to, 0, I) ;
	double ii = ih_interpI(t_to, ind_to, 1, I) ;
	double p[I->length], dt, p_cum = 0;
	double ss = I->t[1] - I->t[0] ;
	int i ;
	double r ;
	double bt ;
	int nstop = n - 1;
	dt = t_from - I->t[ind_from] ;
	
	for (i=ind_from ; i>= ind_to ; i--)
	{
		p_cum += nchoose2 * I->inv_ne[i] * dt ;
		p[i] = exp(-p_cum) ; // probability that coalescent event hasnt occurred   
		r = (double) rand() / (RAND_MAX) ;
		//printf("%8.4f %8.4f %8.4f\n", I->t[i], p[i], r) ;		 	
		if (r > p[i]) // if event has occurred
		{
			p[i] = nchoose2 * I->inv_ne[i] ;
			r = (double) rand() / (RAND_MAX) ;
			r = -log(1 - r*(1 - exp(-dt * p[i]))) / p[i] ; 
			
			if (r < dt)
			{
		 		bt= I->t[i] + dt - r ;
		 		n-- ;
		 		nchoose2 = (n > 0) ? (n * (n-1)) / 2 : 0 ;
		 		p_cum = 0 ;
			}
		}
		
		 
		dt = ss ;
		
		if (n==nstop) break ;
	}
	return(bt) ;
	
}
