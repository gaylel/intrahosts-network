#include "sfs.h"
// functions for calculating (joint) site frequency spectrum 
// represent nodes of phylogeny as a linked list

sfsnode * sfsnode_create(double t, int* desc, int desc_n)
{
	int i ;
	sfsnode *ptr ;
	ptr = calloc(1, sizeof(sfsnode)) ;
	ptr->desc = calloc(desc_n, sizeof(int)) ;
	for (i=0 ; i< desc_n ; i++)
	{
		ptr->desc[i] = desc[i] ; 
	}
	ptr->t = t ;
	ptr->next = (sfsnode *) NULL ;
	return ptr ;
}

void sfsnode_free(sfsnode* x)
{
	free(x->desc) ;
	free(x) ;
	x = NULL ;
}

void sfsnodes_free(sfsnode* cur)
{
	sfsnode* del, *ptr ;
	ptr = cur ;
	while (ptr!=NULL)
	{
		del=ptr->next ;
		sfsnode_free(ptr) ;
		ptr=del;
	}
	
	
}

sfsnode * sfsnode_add(double t, int* desc, int desc_n, sfsnode* cur)
{
	int i ;
	sfsnode *ptr, *head ;
	head = cur ; 
	ptr = sfsnode_create(t, desc, desc_n) ;
	
	if (cur != NULL)
	{
		// existing sfsnode
		while (cur->next!=NULL)
		{
			cur = cur->next ;
		}
		cur->next=ptr ;
		cur = head ;
	}
	else
	{
		cur = ptr ;
	}
	return cur ;
}

int sfsnode_order(int *desc1, int *desc2, int desc_n)
{
	// compare desc1 and desc2 ;
	int i, o = 0 ;
	
	for (i=0 ; i<desc_n ; i++)
	{
		o = desc1[i] - desc2[i] ;
		if (o!=0) break ;
	}
	return o ;
	// if o > 0, desc1 > desc2
	// if o < 0, desc1 < desc2
	// if o==0, desc1 == desc2
	
}


sfsnode * sfsnode_orderedadd(double bl, int *desc, int desc_n, sfsnode* cur)
{
	// add new combination 
	sfsnode *ptr = NULL, *head, *ptr2 ;
	head = cur ; 
	int i , ct = 0, o;
	if (cur == NULL)
	{
		ptr = sfsnode_create(bl, desc, desc_n) ;
		cur = ptr ;
	}
	else
	{
		// iterate 
		o = sfsnode_order(desc, cur->desc, desc_n) ;
		while ((cur->next!=NULL) && (o > 0))
		{
			ptr = cur ;
			cur = cur->next ;
			o = sfsnode_order(desc, cur->desc, desc_n) ;
		}
		
		
		if (o > 0)
		{
			cur->next = sfsnode_create(bl, desc, desc_n) ; 
		
			// desc is current maximum
			// create new
		}
		
		if (o < 0)
		{
			ptr2 = sfsnode_create(bl, desc, desc_n) ; 
			ptr2->next = cur ;
			if (ptr != NULL)
			{
				ptr->next = ptr2 ;
			}
			else
			{
				head = ptr2 ;
			}
		}
		
		if (o == 0)
		{
			cur->t = cur->t + bl ;
		}	
		cur = head ;
	}
	
	return cur ;
}

sfsnode* sfsnode_get(int ind, sfsnode* cur)
{
	sfsnode *ptr = cur;
	int i ;
	for (i=0 ; i<ind ; i++)
	{
		ptr = ptr->next ;
		//cur = ptr ;
	} 
	return ptr ;
}

sfsnode* sfsnode_delete(int ind, sfsnode* cur)
{
	sfsnode *ptr = cur, *del, *head = ptr;
	int i, val ;
	
	// if ind > 0
	for (i=0 ; i<(ind-1) ; i++)
	{
		ptr = cur->next ;
		cur = ptr ;
	} 
	
	// ptr points to entry before ind
	if (ind > 0)
	{
		del=ptr->next ;
		ptr->next = del->next ;
	}
	else
	{
		head = ptr->next ;
		del = ptr ;
	}
	
	free(del->desc) ;
	free(del) ;
	del = NULL ;
	return head ;
}

sfstype* sfs_init(int Nsamp)
{
	sfstype * x;
	x = calloc(1, sizeof(sfstype)) ;
	x->sfs = (sfsnode *) NULL ;
	x->N = Nsamp ;
	x->NNode = 0 ;
	x->sfsp = (sfsnode *) NULL ;
	x->ttot = 0 ;
	return x ;
}

sfstype* sfs_reinit(sfstype* x)
{
	// clear sfs but keep sfsp
	sfsnodes_free(x->sfs) ;
	x->sfs = (sfsnode *) NULL ;
	x->NNode = 0 ;
	return x ;
}

void sfs_free(sfstype* x)
{
	if (x->sfs != NULL) sfsnodes_free(x->sfs) ;
	if (x->sfsp != NULL) sfsnodes_free(x->sfsp) ;
	free(x) ;
} 

sfstype* sfs_addpop(int Nnode, int ind, double t, sfstype* x)
{
	int i ;
	int *desc, desc_n = x->N ; 
	desc = calloc(x->N, sizeof(int)) ;
	// indicate 1 descendant for each node, in host ind.
	for (i=0 ; i<desc_n ; i++)
	{
		desc[i] = 0 ;
	}
	desc[ind] = 1 ;
	
	for (i=0 ; i<Nnode ; i++)
	{
		x->sfs = sfsnode_add(t, desc, desc_n, x->sfs) ;
		x->NNode++ ;
	}
	free(desc) ;
	return x ;
}

sfstype* sfs_coalesce(sfstype* x, double t_c)
{
	// coalesce two of the nodes in x
	// generate random number
	int r, d1[x->N], d2[x->N], i ;
	double t1, t2 ;
	sfsnode* node ;
	
	// first node
	r = (int) ((double)(x->NNode) * rand() / ( RAND_MAX + 1.0));
	//printf("%i\t", r) ;
	node = sfsnode_get(r, x->sfs) ;
	t1 = node->t ;
	for (i=0 ; i<x->N ; i++)
	{
		d1[i] = node->desc[i] ;
	//	printf("%i ", d1[i]) ;
	}
	//printf("\n") ;
	x->sfsp = sfsnode_orderedadd(t1 - t_c, d1, x->N, x->sfsp) ;
	// delete node
	x->sfs = sfsnode_delete(r, x->sfs) ;
	x->NNode-- ;
	
	r = (int) ((double)(x->NNode) * rand() / ( RAND_MAX + 1.0));
	//printf("%i\n", r) ;
	node = sfsnode_get(r, x->sfs) ;
	t2 = node->t ;
	for (i=0 ; i<x->N ; i++)
	{
		d2[i] = node->desc[i] ;
	//	printf("%i ", d2[i]) ;
	}
	//printf("\n") ;
	x->sfsp = sfsnode_orderedadd(t2 - t_c, d2, x->N, x->sfsp) ;
	
	// delete node
	x->sfs = sfsnode_delete(r, x->sfs) ;
	//x->NNode-- ;
	
	// add coalesced node
	for (i=0 ; i<x->N ; i++)
	{
		d1[i] += d2[i] ;
	}
	x->sfs = sfsnode_add(t_c, d1, x->N, x->sfs) ;
	
	return x ;
}

sfstype* sfs_update(int h1, int h2, double t, int k, sfstype* x)
{
	//printf("k=%i", k) ;
	if (k<0)
	{
		// adding nodes
		x = sfs_addpop(-k, h2, t, x) ;
	}
	else
	{
		// coalescent event
		//printf("coalesce %i %8.4f %i %i\n", k, t, h1, h2) ;
		x = sfs_coalesce(x, t) ;
	}
}

double sfs_bltot(sfstype* x)
{
	sfsnode *cur ;
	double ttot = 0;
	cur = x->sfsp ;
	
	while (cur!=NULL)
	{
		ttot+= cur->t ;
		cur = cur->next ;
	}
	return ttot ;
}


double sfs_getbl(sfstype* x, int* desc)
{
	sfsnode *cur ;
	int o = 1;
	double bl = 0 ;
	cur = x->sfsp ;
	
	while (cur!=NULL && o>0)
	{
		o = sfsnode_order(desc, cur->desc, x->N) ;
		if (o == 0) bl = cur->t ;
		cur = cur->next ;
	}
	
	return bl ;
}



void sfs_print(sfsnode* cur, int N)
{
	sfsnode * ptr = cur ;
	int i ;
	
	
	while (ptr!=NULL)
	{
		for (i=0 ; i< N ; i++)
		{
			printf("%i ", ptr->desc[i]) ; 
		}
		printf("%8.4f\n", ptr->t) ; 
				
		ptr = ptr->next ;
		 
	}
	
}
