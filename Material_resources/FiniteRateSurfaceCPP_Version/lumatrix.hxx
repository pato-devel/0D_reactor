/*   
 * lumatrix.hxx
 * ------------
 *    file created:   Febraury   3, 2013
 *    last modified:  Febraury   3, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        LU decomposition functions
 *
 * notable revisions:
 * -----------------
 *    N/A
 *
 *
 */
#include "FRMArray.hxx"


#if !defined(_LU_MATRIX_DOT_H)
#define _LU_MATRIX_DOT_H



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                     PREPROCESSOR DEFINITIONS                           ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###              PUBLIC STRUCTURE (CLASS) DECLARATIONS                     ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC VARIABLE DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC FUNCTION DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/* 
 * name:      lusolve
 *
 * arguments: A  = LHS matrix
 *            b  = RHS matrix
 *            x  = solution matrix
 *
 * return:    0 = no error (success), 1 = array size mismatch
 *
 * purpose:   solve linear system with LU decomposition
 *
 */
int lusolve(Array2D<double> &A, Array1D<double> &b, Array1D<double> &x);



/* 
 * name:      luinverse
 *
 * arguments: A    = original matrix
 *            Ainv = inverse matrix
 *
 * return:    0 = no error (success), 1 = array size mismatch
 *
 * purpose:   invert matrix
 *
 */
int luinverse(Array2D<double> &A, Array2D<double> &Ainv);



/* 
 * name:      ludecomp_pivot
 *
 * arguments: A        = original matrix
 *            lu       = decomposed matrix
 *            indexvec = ordering matrix
 *
 * return:    0 = no error (success), 1 = array size mismatch
 *
 * purpose:   perform LU decomposition with included pivoting
 *
 */ 
int ludecomp_pivot(Array2D<double> &A, Array2D<double> &lu, Array1D<int> &indexvec);



/* 
 * name:      lu_pivot_rhs
 *
 * arguments: indexvec = ordering matrix
 *            rhs      = right hand side matrix to reorder
 *
 * return:    0 = no error (success), 1 = array size mismatch
 *
 * purpose:   re-order the RHS matrix to match the left
 *
 */ 
int lu_pivot_rhs(Array1D<int> &indexvec, Array1D<double> &rhs);



/* 
 * name:      lu_subst
 *
 * arguments: lu  = decomposed matrix
 *            rhs = transformed right hand side
 *            ans = solution matrix
 *
 * return:    0 = no error (success), 1 = array size mismatch
 *
 * purpose:   perform substitution to get solution vector
 *
 */ 
int lu_subst(Array2D<double> &lu, Array1D<double> &rhs, Array1D<double> &ans);


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC OPERATOR DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


#endif   //  _LU_MATRIX_DOT_H
