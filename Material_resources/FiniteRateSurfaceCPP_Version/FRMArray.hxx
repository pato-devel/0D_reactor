/*   
 * FRMArray.hxx
 * ------------
 *    file created:   January   16, 2013
 *    last modified:  March     12, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        define classes to simulate multi-dimensional dynamic arrays
 *
 * notable revisions:
 * -----------------
 *    N/A
 *
 */

#if !defined(_FRM_ARRAY_DOT_H)
#define _FRM_ARRAY_DOT_H

/* global includes */
#include <cstdlib>

#define inline

/* namespace */
using namespace std;


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                     PREPROCESSOR DEFINITIONS                           ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */


//---------------------- DEFINED TYPES ----------------------

typedef unsigned int array_index;

//-------------------- END DEFINED TYPES --------------------


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###              PUBLIC STRUCTURE (CLASS) DECLARATIONS                     ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



/*  ------------------------ <CLASS> Array1D <CLASS> -----------------------------
 *
 *  description:  1D Array class with operators
 *
 *  ------------------------------------------------------------------------------
 */
template <class data_type> class Array1D
{
public:
      //-----constructor-----
      Array1D(array_index row=1);
      //-----copy constructor-----
      Array1D(const Array1D<data_type>& a);
      //-----destructor-----
      ~Array1D();
      
      //-----accessors-----
      inline array_index getRows() const;
      inline size_t      getSizeOf() const;
	
	//-----mutators-----
	inline void initialize(data_type val);
	inline void resize(array_index row);
      
      //-----subscript operators-----
      inline data_type& operator() (array_index row);
      inline data_type  operator() (array_index row) const;
      
      //-----assignment operators-----
	inline Array1D<data_type>& operator= (const Array1D<data_type>& a);
      
protected:
      array_index rows_val;
      data_type *data;

};
/* ---------------------------<END> Array1D <END>--------------------------------- */



/*  ------------------------ <CLASS> Array2D <CLASS> -----------------------------
 *
 *  description:  2D Array class with operators
 *
 *  ------------------------------------------------------------------------------
 */
template <class data_type> class Array2D
{
public:
      //-----constructor-----
      Array2D(array_index row=1, array_index col=1);
      //-----copy constructor-----
      Array2D(const Array2D<data_type>& a);
      //-----destructor-----
      ~Array2D();
      
      //-----accessors-----
      inline array_index getRows() const;
      inline array_index getCols() const;
      inline size_t      getSizeOf() const;
	
	//-----mutators-----
	inline void initialize(data_type val);
	inline void resize(array_index row, array_index col);
      
      //-----subscript operators-----
      inline data_type& operator() (array_index row, array_index col);
      inline data_type  operator() (array_index row, array_index col) const;
      
      //-----assignment operators-----
	inline Array2D<data_type>& operator= (const Array2D<data_type>& a);
      
protected:
      array_index rows_val, cols_val;
      data_type *data;

};
/* ---------------------------<END> Array2D <END>--------------------------------- */



/*  ------------------------ <CLASS> Array3D <CLASS> -----------------------------
 *
 *  description:  3D Array class with operators
 *
 *  ------------------------------------------------------------------------------
 */
template <class data_type> class Array3D
{
public:
      //-----constructor-----
      Array3D(array_index row=1, array_index col=1, array_index sub=1);
      //-----copy constructor-----
      Array3D(const Array3D<data_type>& a);
      //-----destructor-----
      ~Array3D();
      
      //-----accessors-----
      inline array_index getRows() const;
      inline array_index getCols() const;
      inline array_index getSub() const;
      inline size_t      getSizeOf() const;
      
      //-----mutators-----
      inline void setLocation(array_index row = 0, array_index col = 0);
	inline void initialize(data_type val);
	inline void resize(array_index row, array_index col, array_index sub);
      
      //-----subscript operators-----
      inline data_type& operator() (array_index row, array_index col, array_index sub);
      inline data_type  operator() (array_index row, array_index col, array_index sub) const;
      inline data_type& operator() (array_index sub);
      inline data_type  operator() (array_index sub) const;
      
      //-----assignment operators-----
	inline Array3D<data_type>& operator= (const Array3D<data_type>& a);
      
protected:
      array_index rows_val, cols_val, sub_val;
      array_index row_mult, col_mult;
      data_type *data, *location;

};
/* ---------------------------<END> Array3D <END>--------------------------------- */



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
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                   PUBLIC OPERATOR DECLARATIONS                         ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */



#endif   //  _FRM_ARRAY_DOT_H
