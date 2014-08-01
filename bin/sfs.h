#include <stdlib.h>
#include <stdio.h>
#include <time.h>

struct sfsnode {
	int *desc ;
	double t ;
	struct sfsnode * next ;
} ;

typedef struct sfsnode sfsnode;



struct sfstype{
	sfsnode* sfs ;
	int N ;	// number of sampling times
	int NNode ; 
	sfsnode* sfsp ; // record which nodes affected by which branch lengths
	double ttot ; // total branch length in tree
	double *ebt ; // expected branching times
	int ebtl ; //length of ebt ; 
	int ebtp ; // current index into ebt
} ;
typedef struct sfstype sfstype ;

sfsnode * sfsnode_create(double t, int* desc, int desc_n) ;
void sfsnode_free(sfsnode* x) ;
void sfsnodes_free(sfsnode* cur) ;
sfsnode * sfsnode_add(double t, int* desc, int desc_n, sfsnode* cur) ;
sfsnode * sfsnode_add_tofront(double t, int *desc, int desc_n, sfsnode* cur) ;
int sfsnode_order(int *desc1, int *desc2, int desc_n) ;
sfsnode * sfsnode_orderedadd(double bl, int *desc, int desc_n, sfsnode* cur) ;
sfsnode* sfsnode_get(int ind, sfsnode* cur) ;
sfsnode* sfsnode_delete(int ind, sfsnode* cur) ;
sfstype* sfs_update(int h1, int h2, double t, int k, sfstype* x) ;
sfstype* sfs_init(int Nsamp) ;
sfstype* sfs_reinit(sfstype* x) ;
sfstype* sfs_coalesce(sfstype* x, double t_c) ;
sfstype* sfs_addpop(int Nnode, int ind, double t, sfstype* x) ;
void sfs_free(sfstype* x) ;
double sfs_bltot(sfstype* x) ;
double sfs_getbl(sfstype* x, int* desc) ;
void sfs_print(sfsnode* cur, int N) ;
sfstype* sfs_ebtinit(sfstype* x, int* Nseq) ;
sfstype* sfs_ebtadd(sfstype* x, double bt) ;
sfstype* sfs_ebtmean(sfstype* x, int Niter) ;
void sfs_ebtprint(sfstype *x) ;
int sfs_ifexists(int *desc, int desc_n, int**mut, int mut_n) ;
sfstype* sfs_update2(int h1, int h2, double t, int k, sfstype* x, int **mut, int mut_n) ;
sfstype* sfs_coalesce2(sfstype* x, double t_c, int** mut, int mut_n) ;