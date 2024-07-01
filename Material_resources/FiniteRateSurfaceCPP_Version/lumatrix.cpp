/*
 * lumatrix.cpp
 * ------------
 *    file created:   Febraury   3, 2013
 *    last modified:  Febraury   4, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        LU decomposition functions
 *
 * notable revisions:
 * -----------------
 *    N/A
 *
 */
#include <cmath>
#include "lumatrix.hxx"


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                  LOCAL PREPROCESSOR DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                 LOCAL (PRIVATE) VARIABLE DECLARATIONS                  ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###               LOCAL (PRIVATE) FUNCTION DECLARATIONS                    ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                     PUBLIC FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/* -------------------------------- lusolve -------------------------------------- */
/* 
 * solve system of equations: overwrite A and b
 *
 */
int lusolve(Array2D<double> &A, Array1D<double> &b, Array1D<double> &x)
{
	Array1D<int> indexvec;
	int    size;

	/* check that size is consistent */
	size = b.getRows();
	if ( (A.getRows() != size) || (A.getCols() != size) || (x.getRows() != size) )
		return 1;

	/* re-allocate */
	indexvec.resize(size);

	/* LU decomposition */
	ludecomp_pivot(A, A, indexvec);
	/* pivot RHS */
	lu_pivot_rhs(indexvec, b);
	/* substitute */
	lu_subst(A, b, x);
	
	return 0;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- luinverse ----------------------------------- */
/* 
 * compute inverse
 *
 */
int luinverse(Array2D<double> &A, Array2D<double> &Ainv)
{
      Array1D<int> indexvec;
	Array2D<double> lu;
	Array1D<double> rhs, x;
	int c, r, size;
	
	/* check that size is consistent */
	size = A.getRows();
	if ( (A.getCols() != size) || (Ainv.getRows() != size) || (Ainv.getCols() != size) )
		return 1;
	
	/* allocate */
	indexvec.resize(size);
	lu.resize(size,size);
	rhs.resize(size);
	x.resize(size);
      
      /* perform LU decomposition once */
      ludecomp_pivot(A, lu, indexvec);
      
      /* find columns of inverse by successive substitution */
      for (c=0; c<size; c++)  {
            rhs.initialize(0.0);
            rhs(c) = 1.0;
            lu_pivot_rhs(indexvec, rhs);
            lu_subst(lu, rhs, x);
		for (r=0; r<size; r++)
                  Ainv(r,c) = x(r);
      }
	
	return 0;
}
/* ------------------------------------------------------------------------------- */



/* ---------------------------- ludecomp_pivot ----------------------------------- */
/* 
 * perform LU decomposition with included pivoting
 *
 */
int ludecomp_pivot(Array2D<double> &A, Array2D<double> &lu, Array1D<int> &indexvec)
{
      int r, c, k, pindex, itmp, size;
      double pivot,tmp;
	Array1D<double> scale;
	
	/* check that size is consistent */
	size = indexvec.getRows();
	if ( (A.getRows() != size) || (A.getCols() != size) || (lu.getRows() != size) || (lu.getCols() != size) )
		return 1;
      
	scale.resize(size);
      
      /* copy over the original matrix */
      for (r=0; r<size; r++)    {
		for (c=0; c<size; c++)
			lu(r,c) = A(r,c);
	}
      
      /* set index matrix */
      for (c=0; c<size; c++)
            indexvec(c) = c;
      
      /* find largest elements in each row */
      for (r=0; r<size; r++)  {
            scale(r) = 0.0;
            for (c=0; c<size; c++)  {
                  tmp = abs( A(r,c) );
                  if (tmp > scale(r)) scale(r) = tmp;
            }
            scale(r) = 1.0/scale(r);
      }
      
      /* main algorithm */
      for (c=0; c<size; c++)
      {
            /* operate on indices up to not including diagonal */
            for (r=0; r<c; r++)    {
                  for (k=0; k<r; k++)
                        lu(r,c) = lu(r,c) - lu(r,k) * lu(k,c);
            }
            
            pivot = 0.0;
            for (r=c; r<size; r++)    {
                  for (k=0; k<c; k++)
                        lu(r,c) = lu(r,c) - lu(r,k) * lu(k,c);
                  tmp = scale(r) * abs( lu(r,c) );
                  if ( tmp > pivot ) {
                        pivot = tmp;
                        pindex = r;
                  }
            }
            if (pindex != c)  {
                  for (k=0; k<size; k++)     {
                        tmp = lu(pindex,k);
                        lu(pindex,k) = lu(c,k);
                        lu(c,k) = tmp;
                  }
                  itmp = indexvec(c);
                  indexvec(c) = indexvec(pindex);
                  indexvec(pindex) = itmp;
                  scale(pindex) = scale(c);
            }
            for (r=c+1; r<size; r++)
                  lu(r,c) = lu(r,c) / lu(c,c);
      }
      
      return 0;
}
/* ------------------------------------------------------------------------------- */



/* ----------------------------- lu_pivot_rhs ------------------------------------ */
/* 
 * re-order the RHS matrix to match the left
 *
 */
int lu_pivot_rhs(Array1D<int> &indexvec, Array1D<double> &rhs)
{
	Array1D<double> tmp;
	int j, size;

	/* check that size is consistent */
	size = indexvec.getRows();
	if (rhs.getRows() != size)	return 1;

	tmp.resize(size);

	for (j=0; j<size; j++)
		tmp(j) = rhs(indexvec(j));
	for (j=0; j<size; j++)
		rhs(j) = tmp(j);

	return 0;
}
/* ------------------------------------------------------------------------------- */



/* ------------------------------- lu_subst -------------------------------------- */
/* 
 * perform substitution to get solution vector
 *
 */
int lu_subst(Array2D<double> &lu, Array1D<double> &rhs, Array1D<double> &ans)
{
	Array1D<double> ytmp;
	int i, k, size;

	/* check that size is consistent */
	size = rhs.getRows();
	if ( (lu.getRows() != size) || (lu.getCols() != size) || (ans.getRows() != size) )
		return 1;

	ytmp.resize(size);

	/* forward substitution step */
	for (i=0; i<size; i++)  {
		ytmp(i) = rhs(i);
		for (k=0; k<i; k++)
			ytmp(i) = ytmp(i) - lu(i,k) * ytmp(k);
	}
	/* backward substitution step */
	for (i=(size-1); i>=0; i--)  {
		ans(i) = ytmp(i);
		for (k=(i+1); k<size; k++)
			ans(i) = ans(i) - lu(i,k) * ans(k);
		ans(i) = ans(i) / lu(i,i);
	}

	return 0;
}
/* ------------------------------------------------------------------------------- */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                    PRIVATE FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */
