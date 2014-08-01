#include "intrahost2.h"
#include <Rmath.h>


hosttype* ih_hostinit(int* ind, int* Nseq, double *tsamp, int Nsamp, double t_off)
{
	hosttype *h ;
	h = calloc(1, sizeof(hosttype)) ;
	h->no = ind ;
	h->Nsamp = Nsamp ; 
	h->Nseq = Nseq ;
	h->tsamp = tsamp ;
	h->t_off = t_off ; // time to (first) sampling from infection
	h->t_inf = h->tsamp[0] - h->t_off ; // time of (first) infection
	return h ;
}

void ih_hostfree(hosttype *h)
{
	//if (h->no != NULL) free(h->no) ;	// host number
	//if (h->Nseq != NULL) free(h->Nseq) ;	//	number of sequences for each sample
	//if (h->tsamp != NULL) free(h->tsamp) ;	// sampling times
	free(h) ;
}

sfsdatatype* sfsdata_init(int** mut, int *mutc, int mut_n, int lp, int S, double* ts, int *ns)
{
	sfsdatatype* x ; 
	x = calloc(1, sizeof(sfsdatatype)) ;
	x->mut = mut ;
	x->mutc = mutc ;
	x->mut_n = mut_n ;
	x->lp = lp ;
	x->S = S ;
	x->ts = ts ;
	x->ns = ns ;
	return x ;
}

void sfsdata_free(sfsdatatype *x )
{
	free(x->mutc) ;
	free(x->ns) ;
	free(x->ts) ;
	int i ;
	for (i=0 ; i<x->mut_n ; i++)
	{
		free(x->mut[i]) ;
	}
	free(x->mut) ;
	free(x) ;
}


/*
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

void ih_free(Itrajtype* I)
{
	free(I->tc) ;
	free(I->v) ;
	free(I->t) ;
	free(I->i) ;
	free(I->inv_ne) ;
	free(I) ;
}

double* ih_get_cn(hosttype* h, Itrajtype *I)
{
	int ind_from, i ;
	double *cn ;
	cn = calloc(h->Nsamp, sizeof(double)) ;
	for (i=0 ; i<h->Nsamp ; i++)
	{
		ind_from = ih_tlookup(h->tsamp[i] - h->t_inf, I) ;
	
	 	cn[i] = ih_interpI(h->tsamp[i] - h->t_inf, ind_from, 1, I) ;
	}
	return cn ; 
}

int ih_check(double* cn, hosttype *h)
{
	int i , rv=1;
	for (i=0 ; i<h->Nsamp ; i++)
	{
		if (floor(cn[i]) < h->Nseq[i])
		{
			rv = -1 ;
			break ;
		}
	}
	
	return rv ; 
}


Itrajtype * ih_drawI_Baccam(double bet, double delt, double c, double p, double T0, double V0, double ss )
{
	/*******************************************************************************
	
	Draws deterministic trajectory using Euler's forward method
	
	Args:	
	Returns:
			pointer to Itrajtype
	
	********************************************************************************/
	
	double thr = (double ) 1 ;
	double t_i = 0, d_v, d_i, d_tc ;
	int n=0;
	
	Itrajtype* I ;
	I = calloc(1, sizeof(Itrajtype));
	I->length = 1 ;
	I->T0 = T0 ;
	I->V0 = V0 ;
	I->i = calloc(1, sizeof(double)) ; 
	I->i[0] = 0;
	I->tc = calloc(1, sizeof(double)) ; 
	I->tc[0] = T0 ;
	I->v = calloc(1, sizeof(double)) ; 
	I->v[0] = V0;
	I->t = calloc(1, sizeof(double)) ;
	I->t[0] = 0 ;
	
	I->inv_ne = calloc(1, sizeof(double)) ; 
	I->inv_ne[0] = (2 * p * bet * I->tc[0]) / (I->v[0] - 1);
	//I->N = N ;
	I->bet = bet ; 
	I->delt = delt ;
	I->c = c ;
	I->p = p ;
	I->T = 0 ;
	double Tmax = 30.0 ;
	
	while (I->v[n] >= thr && I->T < Tmax)
	{
		d_tc = -(bet * I->tc[n] * I->v[n]) * ss ;
		d_i = -d_tc - (delt * I->i[n] * ss) ;  
		d_v = ((p * I->i[n]) - (c * I->v[n])) * ss ;
		//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", I->t[n], I->v[n], I->tc[n], I->i[n], d_i, d_v, I->T) ;
		
		I->tc = realloc(I->tc,  (I->length + 1) * sizeof(double)) ;
		I->v = realloc(I->v,  (I->length + 1) * sizeof(double)) ;
		I->i = realloc(I->i,  (I->length + 1) * sizeof(double)) ;
		I->t = realloc(I->t,  (I->length + 1) * sizeof(double)) ;
		
		I->inv_ne = realloc(I->inv_ne,  (I->length + 1) * sizeof(double)) ;
		I->length++ ;
		n++ ;
		I->tc[n] = I->tc[n-1] + d_tc ;
		I->tc[n] = (I->tc[n] < 0) ? 0 : I->tc[n] ; 
		I->v[n] = I->v[n-1] + d_v ;
		I->v[n] = (I->v[n] <= 0) ? 0 : I->v[n] ; 
		
		I->i[n] = I->i[n-1] + d_i ;  
		I->i[n] = (I->i[n] <= 0) ? 0 : I->i[n] ; 
		
		I->t[n] = I->t[n-1] + ss ;
		I->inv_ne[n] = (2 * p * bet * I->tc[n]) / (I->v[n] - 1);
		I->T+= ss ; 
		 
	}
	
	return I ;
		
}

Itrajtype * ih_drawI_Baccam2(double bet, double delt, double c, double p, double T0, double V0, double ss )
{
	/*******************************************************************************
	
	Draws deterministic trajectory using Euler's forward method
	
	Args:	
	Returns:
			pointer to Itrajtype
	
	********************************************************************************/
	
	double thr = (double ) 1 ;
	double t_i = 0, d_v, d_i, d_tc ;
	double Tmax = 30.0 ;
	int maxl = (int) ceil(Tmax / ss) ;
	int n=0;
	
	Itrajtype* I ;
	I = calloc(1, sizeof(Itrajtype));
	I->length = 1 ;
	I->T0 = T0 ;
	I->V0 = V0 ;
	
	
	I->i = calloc(maxl, sizeof(double)) ; 
	I->i[0] = 0;
	I->tc = calloc(maxl, sizeof(double)) ; 
	I->tc[0] = T0 ;
	I->v = calloc(maxl, sizeof(double)) ; 
	I->v[0] = V0;
	I->t = calloc(maxl, sizeof(double)) ;
	I->t[0] = 0 ;
	
	I->inv_ne = calloc(maxl, sizeof(double)) ; 
	I->inv_ne[0] = (2 * p * bet * I->tc[0]) / (I->v[0] - 1);
	//I->N = N ;
	I->bet = bet ; 
	I->delt = delt ;
	I->c = c ;
	I->p = p ;
	I->T = 0 ;
	
	
	while (I->v[n] >= thr && I->T < Tmax)
	{
		d_tc = -(bet * I->tc[n] * I->v[n]) * ss ;
		d_i = -d_tc - (delt * I->i[n] * ss) ;  
		d_v = ((p * I->i[n]) - (c * I->v[n])) * ss ;
		//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", I->t[n], I->v[n], I->tc[n], I->i[n], d_i, d_v, I->T) ;
		
		//I->tc = realloc(I->tc,  (I->length + 1) * sizeof(double)) ;
		//I->v = realloc(I->v,  (I->length + 1) * sizeof(double)) ;
		//I->i = realloc(I->i,  (I->length + 1) * sizeof(double)) ;
		//I->t = realloc(I->t,  (I->length + 1) * sizeof(double)) ;
		
		//I->inv_ne = realloc(I->inv_ne,  (I->length + 1) * sizeof(double)) ;
		I->length++ ;
		n++ ;
		I->tc[n] = I->tc[n-1] + d_tc ;
		I->tc[n] = (I->tc[n] < 0) ? 0 : I->tc[n] ; 
		I->v[n] = I->v[n-1] + d_v ;
		I->v[n] = (I->v[n] <= 0) ? 0 : I->v[n] ; 
		
		I->i[n] = I->i[n-1] + d_i ;  
		I->i[n] = (I->i[n] <= 0) ? 0 : I->i[n] ; 
		
		I->t[n] = I->t[n-1] + ss ;
		I->inv_ne[n] = (2 * p * bet * I->tc[n]) / (I->v[n] - 1);
		I->T+= ss ; 
		 
	}
	
	return I ;
		
}


ivectype* ih_Ilookup(int l, Itrajtype* I)
{
/*******************************************************************************
	
	Gets index of deterministic trajectory for value of I
	
	Args:	l number of sequences
			I Itraj
			
	Returns:
			index of I 
	
	********************************************************************************/
	int n, nfound=0, *ind ;
	ivectype* v = calloc(1, sizeof(ivectype)) ;
	for (n=1 ; n< (I->length - 1) ; n++)
	{
		if ((((I->v[n]) > (double) l) && ((I->v[n-1]) <= (double) l)) || (((I->v[n+1]) <= (double) l) && ((I->v[n]) > (double) l)))
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
	double diff, frac_t, rv=-1 ;
	
	if (((t_i + 1) < I->length) && (t_i>=0))
	{
	frac_t = I->t[t_i+1] - I->t[t_i] ;
	frac_t = (t - I->t[t_i]) / frac_t ;   
	switch(isI)
	{
		case 0: // tc
			diff = I->tc[t_i + 1] - I->tc[t_i] ;
			rv = I->tc[t_i] + frac_t * diff ; 	 
		break ;
		case 1: // v
			diff = I->v[t_i + 1] - I->v[t_i] ;
			rv = I->v[t_i] + frac_t * diff ;
		case 2: // i
			diff = I->i[t_i + 1] - I->i[t_i] ;
			rv = I->i[t_i] + frac_t * diff ;
		break ;
	}
	}
	return rv ;
}

/*
void ih_fit_pair(Itrajtype *I1, Itrajtype *I2, int n1, int n2, double t1, double t2)
{
	// scale times
	// t2 and t1 on traj. time	

	// draw coalescent for recent host
	int Nav = 1 ;
	int i, nlin2=n2 ;
	double *bt2, *bt1 ; 
	
	bt2 = ih_draw_ebt(I2, n2, Nav, t2, t1) ;
	
	// see how many lineages left at the start;
	for (i=(n2 - 2 ) ; i>= 0 ; i--)
	{
		if (bt2[i] <= 0)
		{
			break ;
		}
		nlin2 -- ;
	}
	
	bt1 = ih_draw_ebt(I1, n1 + nlin2, Nav, t1, 0) ;
	// 
	
	double **jsfs = calloc(n1, sizeof(double*)) ;
	for (i=0 ; i< n1 ; i++)
	{
		jsfs[i] = calloc(n2, sizeof(double)) ;
	}
	
	// draw tree
	
	 
	// 
}*/

bttype * ih_btmxinit()
{
	bttype *bt ;
	bt = calloc(1, sizeof(bttype)) ;
	/*int *nlin ;	// number of lineages
	int *h1 ; // host index
	int *h2 ;
	double *t ;	// branching time 
	int *k ; // key 1=coalescent, -n=adding n lineages
	*/
	bt->N = 0;
	return bt ;
}

bttype * ih_btmxinit2(int* Nseq, int Nsamp)
{
	bttype *bt ;
	bt = calloc(1, sizeof(bttype)) ;
	/*int *nlin ;	// number of lineages
	int *h1 ; // host index
	int *h2 ;
	double *t ;	// branching time 
	int *k ; // key 1=coalescent, -n=adding n lineages
	*/
	
	int i, nlin=0 ;
	for (i = 0 ; i< Nsamp ; i++)
	{
		nlin += Nseq[i] ;
	}
	bt->nlin = calloc(nlin, sizeof(int)) ;
	bt->h1 = calloc(nlin, sizeof(int)) ;
	bt->h2 = calloc(nlin, sizeof(int)) ;
	bt->k = calloc(nlin, sizeof(int)) ;
	bt->t = calloc(nlin, sizeof(double)) ;
	bt->N = 0;
	return bt ;
}


bttype * ih_btmxreset(bttype *x)
{
	free(x->nlin) ;
	free(x->h1) ;
	free(x->h2) ;
	free(x->t) ;
	free(x->k) ;
	x->N = 0 ;
	
	x->nlin = (int*) NULL ;
	x->h1 = (int*) NULL ;
	x->h2 = (int*) NULL ;
	x->k = (int*) NULL ;
	x->t = (double*) NULL ;
	return x ;	
	
}

void ih_btfree(bttype *x) 
{
	if (x->nlin != NULL) free(x->nlin) ;
	if (x->h1 != NULL) free(x->h1) ;
	if (x->h2 != NULL) free(x->h2) ;
	if (x->t != NULL) free(x->t) ;
	if (x->k != NULL) free(x->k) ;
	free(x) ;
}



bttype * ih_btmxupdate(bttype* bt, int nlin, int h1, int h2, double t, int k)
{
	if (bt->N == 0)
	{
		bt->nlin = calloc(1, sizeof(int)) ;
		bt->h1 = calloc(1, sizeof(int)) ;
		bt->h2 = calloc(1, sizeof(int)) ;
		bt->k = calloc(1, sizeof(int)) ;
		bt->t = calloc(1, sizeof(double)) ;
	}
	else
	{
		bt->nlin = realloc(bt->nlin, (bt->N+1) *sizeof(int)) ;
		bt->h1 = realloc(bt->h1, (bt->N+1) * sizeof(int)) ;
		bt->h2 = realloc(bt->h2, (bt->N+1) * sizeof(int)) ;
		bt->k = realloc(bt->k, (bt->N+1) * sizeof(int)) ;
		bt->t = realloc(bt->t, (bt->N+1) * sizeof(double)) ;
		
	}
	bt->nlin[bt->N] = nlin ;
	bt->h1[bt->N] = h1 ;
	bt->h2[bt->N] = h2 ;
	bt->k[bt->N] = k ;
	bt->t[bt->N] = t;
	//printf("%i %i %i %8.4f %i\n", nlin, h1, h2, t, k) ; 
	bt->N++ ;
	return bt ;
}

bttype * ih_btmxupdate2(bttype* bt, int nlin, int h1, int h2, double t, int k)
{
	bt->nlin[bt->N] = nlin ;
	bt->h1[bt->N] = h1 ;
	bt->h2[bt->N] = h2 ;
	bt->k[bt->N] = k ;
	bt->t[bt->N] = t;
	//printf("%i %i %i %8.4f %i\n", nlin, h1, h2, t, k) ; 
	bt->N++ ;
	return bt ;
}

double* ih_get_ebt(sfsdatatype* sdt, Itrajtype *I, hosttype *h, int Niter)
{
	// get expected branching times of coalescent
	sfstype *x ;
	x = ih_draw_allbt(Niter, h, I) ;
	x = sfs_ebtmean(x, Niter) ;
	//x = ih_draw_coaltree2(I, h, Niter, sdt) ;
	// copy branching times
	double *ebt ;
	ebt = calloc(x->ebtl, sizeof(double)) ;
	int i ;
	for (i=0 ; i<x->ebtl ; i++)
	{
		ebt[i] = x->ebt[i] ; 
	}
	sfs_free(x) ;
	return ebt ;
}


double* ih_psfs(sfsdatatype* sdt, Itrajtype *I, hosttype *h, int Niter)
{
	sfstype * x ;
	double *ll ;
	//x = ih_draw_coaltree(I, h, Niter) ;
	//ll = sfs_calcll(x, sdt, Niter) ;	
	//sfs_print(x->sfsp, x->N) ;
	
	x = ih_draw_coaltree2(I, h, Niter, sdt) ;
	ll = sfs_calcll2(x, sdt, Niter) ;	
	//sfs_print(x->sfsp, x->N) ;
	
	
	sfs_free(x) ;

	return ll ;
}

double* ih_psfs2(sfsdatatype* d, Itrajtype *I, hosttype *h, int Niter)
{
	
	double *ll ;
	double ttot=0 ;
	int i, j ;
	
	sfstype *x ;
	x = ih_draw_allbt(Niter, h, I) ;
	
	//for (i=0 ; i<x->ebtl ; i++)
	//{
	//	ttot += x->ebt[i] ; 
	//}
	//ttot = ttot / Niter ; 
	x = sfs_ebtmean(x, Niter) ;
	//sfs_ebtprint(x) ;
	
	double sumkET = ih_calculateess1(x->ebt, (x->ebtl + 1)) ;
	
	
	double p ;
	ll = calloc(2, sizeof(double)) ;
	ll[0] = 0.0 ;
	ll[1] = sumkET  ;
	for (i=0 ; i<d->mut_n ; i++)
	{
		p = ih_calculateess2(x->ebt, x->ebtl + 1, d->mut[i][0], sumkET) ;
		//printf("b=%i, p=%8.4f, ", d->mut[i][0], p) ;
		ll[0]+= d->mutc[i]*(log(p)) ;
	}
	//printf("\n") ;
	//printf("ttot = %8.4f\n", ttot) ;
	
	// add multinomial term
	ll[0]+= lgammafn(d->mut_n + 1) ;
	for (i=0 ; i<d->mut_n ; i++)
	{
		ll[0]-=  lgammafn(d->mutc[i] + 1) ;
	}
	
	sfs_free(x) ;
	return ll ;
}



sfstype* ih_draw_coaltree2(Itrajtype *I, hosttype* h, int Niter, sfsdatatype *sdt)
{
	double t_from ;  
	double t_to ; 
	double next_bt, old_bt ; 
	int i, nlin, k ; 
	bttype* bt ; 
	sfstype * x ;
	int j ;
	x = sfs_init(h->Nsamp) ;
	sfsnode *node ; 
	int d1[x->N], btn ;
	double t1 ;
	
	bt = ih_btmxinit() ;
	x = sfs_ebtinit(x, h->Nseq) ;
	btn = x->ebtl + h->Nsamp ; 
	bt->nlin = calloc(btn, sizeof(int)) ;
	bt->h1 = calloc(btn, sizeof(int)) ;
	bt->h2 = calloc(btn, sizeof(int)) ;
	bt->k = calloc(btn, sizeof(int)) ;
	bt->t = calloc(btn, sizeof(double)) ;
			
	for (j=0; j< Niter; j++)
	{
		nlin = 0 ;
		old_bt = h->tsamp[h->Nsamp-1] - h->t_inf;
		x->ebtp = x->ebtl - 1 ;
		
		for (i=h->Nsamp-1 ; i>=0 ; i--)
		{
			t_from = h->tsamp[i] - h->t_inf ;
			t_to = (i > 0) ? (h->tsamp[i - 1] - h->t_inf) : 0 ;
			next_bt = t_from ;
			nlin += h->Nseq[i] ;
			//printf("%8.4f\t%8.4f\t%8.4f\n", t_from, t_to, h->t_inf) ;
			bt = ih_btmxupdate2(bt, nlin, h->no[i], h->no[i], next_bt + h->t_inf, -h->Nseq[i]) ;
			while ((next_bt > t_to))
			{
				next_bt = ih_draw_nextbt(I, nlin, next_bt, t_to) ;
				//if (next_bt==-1) i=0 ;
				
				if (next_bt > t_to)
				{
					x = sfs_ebtadd(x, old_bt - next_bt) ;
					old_bt = next_bt ;

					//printf("%i\t%8.4f\n", nlin, next_bt) ;
					nlin-- ;
					bt = ih_btmxupdate2(bt, nlin, h->no[i], h->no[i], next_bt + h->t_inf, 1) ;
					
				}
			
			}
			//
			
	       if (nlin == 2 && i==0) break ;
		}
	
	
	
		for (i=0 ; i<bt->N ; i++)
		{
			//printf("%i of %i :  %8.4f  %i\n", i, bt->N, bt->t[i], bt->k[i]) ;
			x = sfs_update2(bt->h1[i], bt->h2[i], bt->t[i], bt->k[i], x, sdt->mut, sdt->mut_n) ;
		}
		// check to see if mrca is after time of infection
		
		
		// update x->sfsp
		node = x->sfs ; 
		for (i=0 ; i<x->NNode ; i++)
		{
			//node = sfsnode_get(i, x->sfs) ;
			
			for (k=0 ; k<x->N ; k++)
			{
				d1[k] = node->desc[k] ;
				//	printf("%i ", d1[j]) ;
			}
			t1 = node->t ;
			//printf("%8.4f\n", t1) ;
			
			x->ttot += t1 - h->t_inf ;	
			if (sfs_ifexists(d1, x->N, sdt->mut, sdt->mut_n) == 1)
				x->sfsp = sfsnode_orderedadd(t1 - h->t_inf, d1, x->N, x->sfsp) ;
			node = node->next ; 
		}
	

		
		
		
		
		//sfs_print(x->sfsp, x->N) ;
		x = sfs_reinit(x) ;	
		
		bt->N = 0 ;
		//bt = ih_btmxreset(bt) ;
	}

	
		
		
	


	//sfs_calcll(x, sfsdatatype *d, int Niter)
	ih_btfree(bt) ;
	bt = NULL ;
	x = sfs_ebtmean(x, Niter) ;
	//sfs_ebtprint(x) ;
	return x ;
}


sfstype* ih_draw_allbt(int Niter, hosttype* h, Itrajtype* I)
{
	double t_from ;  
	double t_to ; 
	double next_bt, old_bt ; 
	int i, nlin, k ; 
	sfstype * x ;
	int j ;
	x = sfs_init(h->Nsamp) ;
	sfsnode *node ; 
	int d1[x->N], btn ;
	double t1 ;
	
	x = sfs_ebtinit(x, h->Nseq) ;
			
	for (j=0; j< Niter; j++)
	{
		nlin = 0 ;
		old_bt = h->tsamp[h->Nsamp-1] - h->t_inf;
		x->ebtp = x->ebtl - 1 ;
		
		for (i=h->Nsamp-1 ; i>=0 ; i--)
		{
			t_from = h->tsamp[i] - h->t_inf ;
			t_to = (i > 0) ? (h->tsamp[i - 1] - h->t_inf) : 0 ;
			next_bt = t_from ;
			nlin += h->Nseq[i] ;
			while ((next_bt > t_to))
			{
				next_bt = ih_draw_nextbt(I, nlin, next_bt, t_to) ;
				//if (next_bt==-1) i=0 ;
				
				if (next_bt > t_to)
				{
					x = sfs_ebtadd(x, old_bt - next_bt) ;
					old_bt = next_bt ;

					//printf("%i\t%8.4f\n", nlin, next_bt) ;
					nlin-- ;
				}
			
			}
			//
			
	       if (nlin == 2 && i==0) break ;
		}
	}

	
		
		
	


	//x = sfs_ebtmean(x, Niter) ;
	//sfs_ebtprint(x) ;
	return x ;
}

sfstype* ih_draw_coaltree(Itrajtype *I, hosttype* h, int Niter)
{
	double t_from ;  
	double t_to ; 
	double next_bt, old_bt ; 
	int i, nlin, k ; 
	bttype* bt ; 
	sfstype * x ;
	int j ;
	x = sfs_init(h->Nsamp) ;
	sfsnode *node ; 
	int d1[x->N], btn ;
	double t1 ;
	
	bt = ih_btmxinit() ;
	x = sfs_ebtinit(x, h->Nseq) ;
	btn = x->ebtl + h->Nsamp ; 
	bt->nlin = calloc(btn, sizeof(int)) ;
	bt->h1 = calloc(btn, sizeof(int)) ;
	bt->h2 = calloc(btn, sizeof(int)) ;
	bt->k = calloc(btn, sizeof(int)) ;
	bt->t = calloc(btn, sizeof(double)) ;
			
	for (j=0; j< Niter; j++)
	{
		nlin = 0 ;
		old_bt = h->tsamp[h->Nsamp-1] - h->t_inf;
		x->ebtp = x->ebtl - 1 ;
		
		for (i=h->Nsamp-1 ; i>=0 ; i--)
		{
			t_from = h->tsamp[i] - h->t_inf ;
			t_to = (i > 0) ? (h->tsamp[i - 1] - h->t_inf) : 0 ;
			next_bt = t_from ;
			nlin += h->Nseq[i] ;
			//printf("%8.4f\t%8.4f\t%8.4f\n", t_from, t_to, h->t_inf) ;
			bt = ih_btmxupdate2(bt, nlin, h->no[i], h->no[i], next_bt + h->t_inf, -h->Nseq[i]) ;
			while ((next_bt > t_to))
			{
				next_bt = ih_draw_nextbt(I, nlin, next_bt, t_to) ;
				//if (next_bt==-1) i=0 ;
				
				if (next_bt > t_to)
				{
					x = sfs_ebtadd(x, old_bt - next_bt) ;
					old_bt = next_bt ;

					//printf("%i\t%8.4f\n", nlin, next_bt) ;
					nlin-- ;
					bt = ih_btmxupdate2(bt, nlin, h->no[i], h->no[i], next_bt + h->t_inf, 1) ;
					
				}
			
			}
			//
			
	       if (nlin == 2 && i==0) break ;
		}
	
	
	
		for (i=0 ; i<bt->N ; i++)
		{
			//printf("%i of %i :  %8.4f  %i\n", i, bt->N, bt->t[i], bt->k[i]) ;
			x = sfs_update(bt->h1[i], bt->h2[i], bt->t[i], bt->k[i], x) ;
		}
		// check to see if mrca is after time of infection
		
		
		
		node = x->sfs ; 
		for (i=0 ; i<x->NNode ; i++)
		{
			//node = sfsnode_get(i, x->sfs) ;
			
			for (k=0 ; k<x->N ; k++)
			{
				d1[k] = node->desc[k] ;
				//	printf("%i ", d1[j]) ;
			}
			t1 = node->t ;
			//printf("%8.4f\n", t1) ;
			x->sfsp = sfsnode_orderedadd(t1 - h->t_inf, d1, x->N, x->sfsp) ;
			node = node->next ; 
		}
	

		
		
		
		
		//sfs_print(x->sfsp, x->N) ;
		x = sfs_reinit(x) ;	
		
		bt->N = 0 ;
		//bt = ih_btmxreset(bt) ;
	}

	
		
		
	


	//sfs_calcll(x, sfsdatatype *d, int Niter)
	ih_btfree(bt) ;
	bt = NULL ;
	x = sfs_ebtmean(x, Niter) ;
	//sfs_ebtprint(x) ;
	return x ;
}

/*void ih_draw_coalevent(bttype* bt, int n)
{
	
}
*/


/*
double ih_draw_nextct(Itrajtype* I, int NHosts, int** nlin, double *nlin_t, double t_now, double t_min)
{
	double ct, last = t_now ; 


	ih_draw_nextbt(I, int n, double t_now, t_min)
}
*/
//double* ih_draw_sfs(Itrajtype *I, int Nav, double )


/*double* ih_draw_btimes(Itrajtype*, int* ns, double *ts, int NHosts, double mt)
{
	// assume pair for now
	int i ;
	int *nlin = calloc(NHosts, sizeof(int)) ; // number of lineages in each host
	
	// get times of infection
	
	ind_from = ih_tlookup(ts[0], I) ;
	
}*/


double* ih_draw_ebt(Itrajtype* I, int n, int N, double t_from, double t_to)
{

	double *bt =calloc(N, sizeof(double)) , *ebt = calloc(n-1, sizeof(double));
	int ni = n, i, j;
	double last = t_from ; 
	int bt_exit = 1 ;
	double ij ;
	// check to see if params are valid
	int ind_from ;
	if (t_from < I->T )
	{
		ind_from = ih_tlookup(t_from, I) ;
		ij = ih_interpI(t_from, ind_from, 1, I) ;
		if (floor(ij) >= n)
		{
			bt_exit = 0 ;
		}
	}
	
	for (j=n-2 ; j>=0 ; j--)
	{
		ebt[j] = 0 ;
	}
	
	if (bt_exit == 0)
	{
	for (j=n-2 ; j>=0 ; j--)
	{
		//ebt[j] = 0 ;
	for (i=0 ; i<N ; i++)
	{
		bt[i] = ih_draw_nextbt(I, j + 2, last, t_to) ;
		ebt[j] += bt[i] ;
		//printf("%8.4f\n", bt[i]) ;
	}
		ebt[j] = ebt[j] / N ;
		ebt[j] = last - ebt[j] ;
		if (ebt[j] < 0 || ((last - ebt[j]) < 0))
		{
			ebt[j] = 0 ;
			break ;
		}
		last = last - ebt[j] ; 
		 //printf("%8.4f\t%8.4f\n", ebt[j], last) ;
	}
	}
	else
	{
		ebt[0] = -1 ;
	}
	free(bt) ; 
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


double ih_calculateess1(double *ebt, int n)
{
	// calculate sumkET
	int L = n-1, i ;
	double sumkET = 0 ;
	for (i=0 ; i<L ; i++)
	{
		sumkET += (i+2) * ebt[i] ; 
	}
	return sumkET ;	
}

double ih_calculateess2(double *ebt, int n, int i, double sumkET)
{
	// calculate expected site spectrum
	int L = n -1, k;
	
	
	double p = 0 ;
	for (k=2 ; k<= (n - i + 1) ; k++)
	{
		p += k * (k-1) *choose(n-k, i-1) * ebt[k-2] ;
	}
	p = p * beta(n-i,i) / sumkET ;
	return p;
}

double ih_nmutdist(double *ebt, int n)
{
	// approximate the distribution over the number of mutations as poisson with mean sumKET * mutation rate
	// return poisson mean.
	
	double sumkET = 0 ;
	int i, L= n - 1;
	for (i=0 ; i<L ; i++)
	{
		sumkET += (i+2) * ebt[i] ; 
	}
	return sumkET ; 
}

double ih_draw_nextbt(Itrajtype* I, int n, double t_from, double t_to)
{
	// draw expected coalescence time
	int ind_from, ind_to ;
	int nchoose2 = (n * (n-1)) / 2 ;
	ind_from = ih_tlookup(t_from, I) ;
	double bt = -1;
	//double sj = ih_interpI(t_from, ind_from, 0, I) ;
	//double ij = ih_interpI(t_from, ind_from, 1, I) ;
	ind_to = ih_tlookup(t_to, I) ;
	//double si = ih_interpI(t_to, ind_to, 0, I) ;
	//double ii = ih_interpI(t_to, ind_to, 1, I) ;
	double p[I->length], dt, p_cum = 0;
	double ss = I->t[1] - I->t[0] ;
	int i ;
	double r ;
	
	int nstop = n - 1;
	dt = t_from - I->t[ind_from] ;
	if (n>1)
	{
	for (i=ind_from ; i>= ind_to ; i--)
	{
		p_cum += nchoose2 * I->inv_ne[i] * dt ;
		p[i] = exp(-p_cum) ; // probability that coalescent event hasnt occurred   
		r = (double) rand() / (RAND_MAX) ;
		//printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", I->t[i], p[i], r, dt, p_cum, (I->v[i]*I->v[i]) * I->inv_ne[i]) ;		 	
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
		 		//printf("%8.4f %i\n", bt, n, i) ;
			}
		}
		
		 
		dt = ss ;
		
		if (n==nstop) break ;
	}
	}
	return(bt) ;
	
}

double* sfs_calcll(sfstype* x, sfsdatatype *d, int Niter)
{
	double ttot = sfs_bltot(x) , lttot = log(ttot);
	double p ;
	double *ll ;
	ll = calloc(2, sizeof(double)) ;
	ll[0] = 0.0 ;
	ll[1] = ttot / Niter ;
	int i, j;
	for (i=0 ; i<d->mut_n ; i++)
	{
		p = sfs_getbl(x, d->mut[i]) ;
		/*for (j=0 ; j<d->lp ; j++)
		{
			printf("%i ", d->mut[i][j]) ;
		}
		printf("%8.4f\n", p) ;*/
		ll[0]+= d->mutc[i]*(log(p) - lttot) ;
	}
	
	// add multinomial term
	ll[0]+= lgammafn(d->mut_n + 1) ;
	for (i=0 ; i<d->mut_n ; i++)
	{
		ll[0]-=  lgammafn(d->mutc[i] + 1) ;
	}
	
	return ll ;
}

double* sfs_calcll2(sfstype* x, sfsdatatype *d, int Niter)
{
	double ttot = x->ttot , lttot = log(ttot);
	double p ;
	double *ll ;
	ll = calloc(2, sizeof(double)) ;
	ll[0] = 0.0 ;
	ll[1] = ttot / Niter ;
	int i, j;
	for (i=0 ; i<d->mut_n ; i++)
	{
		p = sfs_getbl(x, d->mut[i]) ;
		/*for (j=0 ; j<d->lp ; j++)
		{
			printf("%i ", d->mut[i][j]) ;
		}
		printf("%8.4f\n", p) ;*/
		ll[0]+= d->mutc[i]*(log(p) - lttot) ;
	}
	
	// add multinomial term
	ll[0]+= lgammafn(d->mut_n + 1) ;
	for (i=0 ; i<d->mut_n ; i++)
	{
		ll[0]-=  lgammafn(d->mutc[i] + 1) ;
	}
	
	return ll ;
}


