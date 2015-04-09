#include "main.h"
#include "kd_op.h"

void print_kd(struct perm_pt *p_pts, struct kdnode *p)
{
	if( (p->loson == NULL) || (p->hison == NULL) )
	{
		printf("leaf %d:%d %d\n", p->lopt, p->hipt, p->bucket);
		return;
	}
	else
	{
		printf("%d:%d\n", p->cutdim, p->cutval);
		if( p != NULL  )
		{
			print_kd(p_pts, p->loson);
			print_kd(p_pts, p->hison);
		}
		else return;
	}
}

int find_pred_blk(struct kdnode *p, int x, int y)
{
	int res;

	if( p->bucket == 1 ) 
	{
		return(p->lopt);
	}
	else
	{
		if( p->cutdim == 1 )
		{
			if( x <= p->cutval ) res = find_pred_blk(p->loson, x, y);
			else res = find_pred_blk(p->hison, x, y);
		}
		else 
		{
			if( y <= p->cutval ) res = find_pred_blk(p->loson, x, y);
			else res = find_pred_blk(p->hison, x, y);
		}
		return(res);
	}
}

int find_successor(struct kdnode *p, int xval, int yval)
{
	int res; 

	if( p->bucket == 1) 
	{
		return(p->hipt);
	}
	else
	{
		if( p->cutdim == 1 )
		{
			if( xval < p->cutval ) res = find_successor(p->loson, xval, yval);
			else res = find_successor(p->hison, xval, yval);
		}
		else 
		{
			if( yval < p->cutval ) res = find_successor(p->loson, xval, yval);
			else res = find_successor(p->hison, xval, yval);
		}
		return(res);
	}
}
