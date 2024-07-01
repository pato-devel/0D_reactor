/*
 * FRMArray.cpp
 * ------------
 *    file created:   January   16, 2013
 *    last modified:  March     12, 2013
 *    author:         Matthew MacLean
 *                    maclean@cubrc.org
 *    purpose:        implementation of template class functions for Multidimensional Arrays
 *
 * notable revisions:
 * -----------------
 *    N/A
 *
 *
 */
#include "FRMArray.hxx"
#include <string>


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


//*********************************************
//  Array1D member functions
//*********************************************
template <class data_type> Array1D<data_type>::Array1D(array_index row)
{
      data = NULL;
	resize(row);
}


template <class data_type> Array1D<data_type>::Array1D(const Array1D<data_type>& a)
{
      rows_val = a.getRows();
      data = new data_type[rows_val];
      
      for (array_index i=0; i<rows_val; i++)
		(*this)(i) = a(i);
}


template <class data_type> Array1D<data_type>::~Array1D()
{
      delete [] data;
}


template <class data_type> inline array_index Array1D<data_type>::getRows() const
{     return rows_val;   }


template <class data_type> inline size_t Array1D<data_type>::getSizeOf() const
{     return ( rows_val*sizeof(data_type) + sizeof(array_index) );    }


template <class data_type> inline void Array1D<data_type>::initialize(data_type val)
{
      for (array_index i=0; i<rows_val; i++)
		(*this)(i) = val;
}


template <class data_type> inline void Array1D<data_type>::resize(array_index row)
{
      if (data != NULL) delete [] data;
	rows_val = row;
	if (rows_val <= 0) rows_val = 1;
      data = new data_type[rows_val];
}


template <class data_type> inline data_type& Array1D<data_type>::operator() (array_index row)
{     return *(data + row);  }


template <class data_type> inline data_type  Array1D<data_type>::operator() (array_index row) const
{     return *(data + row);  }


template <class data_type> inline Array1D<data_type>& Array1D<data_type>::operator= (const Array1D<data_type>& a)
{
	this->resize(a.getRows());

      for (array_index i=0; i<a.getRows(); i++)
		(*this)(i) = a(i);

	return *this;
}


//*********************************************
//  Array2D member functions
//*********************************************
template <class data_type> Array2D<data_type>::Array2D(array_index row, array_index col)
{
      data = NULL;
	resize(row, col);
}


template <class data_type> Array2D<data_type>::Array2D(const Array2D<data_type>& a)
{
      rows_val = a.getRows();
      cols_val = a.getCols();
      data = new data_type[rows_val * cols_val];
      
      for (array_index i=0; i<rows_val; i++)
      {
            for (array_index j=0; j<cols_val; j++)
                  (*this)(i,j) = a(i,j);
      }
}


template <class data_type> Array2D<data_type>::~Array2D()
{
      delete [] data;
}


template <class data_type> inline array_index Array2D<data_type>::getRows() const
{     return rows_val;   }


template <class data_type> inline array_index Array2D<data_type>::getCols() const
{     return cols_val;   }


template <class data_type> inline size_t Array2D<data_type>::getSizeOf() const
{     return ( rows_val*cols_val*sizeof(data_type) + 2*sizeof(array_index) );    }


template <class data_type> inline void Array2D<data_type>::initialize(data_type val)
{
      for (array_index i=0; i<rows_val; i++)
      {
            for (array_index j=0; j<cols_val; j++)
                  (*this)(i,j) = val;
      }
}


template <class data_type> inline void Array2D<data_type>::resize(array_index row, array_index col)
{
      if (data != NULL) delete [] data;
	rows_val = row;
      cols_val = col;
	if (rows_val <= 0) rows_val = 1;
	if (cols_val <= 0) cols_val = 1;
      data = new data_type[rows_val * cols_val];
}


template <class data_type> inline data_type& Array2D<data_type>::operator() (array_index row, array_index col)
{     return *(data+row*cols_val + col);  }


template <class data_type> inline data_type  Array2D<data_type>::operator() (array_index row, array_index col) const
{     return *(data+row*cols_val + col);  }


template <class data_type> inline Array2D<data_type>& Array2D<data_type>::operator= (const Array2D<data_type>& a)
{
	this->resize(a.getRows(), a.getCols());

      for (array_index i=0; i<a.getRows(); i++)
      {
            for (array_index j=0; j<a.getCols(); j++)
                  (*this)(i,j) = a(i,j);
      }

	return *this;
}


//*********************************************
//  Array3D member functions
//*********************************************

template <class data_type> Array3D<data_type>::Array3D(array_index row, array_index col, array_index sub)
{
      data = NULL;
	resize(row, col, sub);
}


template <class data_type> Array3D<data_type>::Array3D(const Array3D<data_type>& a)
{
      rows_val = a.getRows();
      cols_val = a.getCols();
      sub_val  = a.getSub();
      row_mult = sub_val*cols_val;
      col_mult = sub_val;
      data = new data_type[rows_val * cols_val * sub_val];
      location = data;
      
      for (array_index i=0; i<rows_val; i++)
      {
            for (array_index j=0; j<cols_val; j++)
            {
                  for (array_index m=0; m<sub_val; m++)
                        (*this)(i,j,m) = a(i,j,m);
            }
      }
}


template <class data_type> Array3D<data_type>::~Array3D()
{
      location = NULL;
      delete [] data;
}


template <class data_type> inline void Array3D<data_type>::setLocation(array_index row, array_index col)
{
      location = data+row*row_mult + col*col_mult;
}


template <class data_type> inline array_index Array3D<data_type>::getRows() const
{     return rows_val;   }


template <class data_type> inline array_index Array3D<data_type>::getCols() const
{     return cols_val;   }


template <class data_type> inline array_index Array3D<data_type>::getSub() const
{     return sub_val;   }


template <class data_type> inline size_t Array3D<data_type>::getSizeOf() const
{     return ( sub_val*rows_val*cols_val*sizeof(data_type) + 5*sizeof(array_index) );    }


template <class data_type> inline void Array3D<data_type>::initialize(data_type val)
{
      for (array_index i=0; i<rows_val; i++)
      {
            for (array_index j=0; j<cols_val; j++)
            {
                  for (array_index m=0; m<sub_val; m++)
                        (*this)(i,j,m) = val;
            }
      }
}


template <class data_type> inline void Array3D<data_type>::resize(array_index row, array_index col, array_index sub)
{
      if (data != NULL) delete [] data;
	rows_val = row;
      cols_val = col;
      sub_val  = sub;
	if (rows_val <= 0) rows_val = 1;
	if (cols_val <= 0) cols_val = 1;
	if (sub_val <= 0)  sub_val  = 1;
      row_mult = sub_val*cols_val;
      col_mult = sub_val;
      data = new data_type[rows_val * cols_val * sub_val];
      location = data;
}


template <class data_type> inline data_type& Array3D<data_type>::operator() (array_index row, array_index col, array_index sub)
{     return *(data+row*row_mult + col*col_mult + sub);  }


template <class data_type> inline data_type  Array3D<data_type>::operator() (array_index row, array_index col, array_index sub) const
{     return *(data+row*row_mult + col*col_mult + sub);  }


template <class data_type> inline data_type& Array3D<data_type>::operator() (array_index sub)
{     return *(location + sub);  }


template <class data_type> inline data_type  Array3D<data_type>::operator() (array_index sub) const
{     return *(location + sub);  }


template <class data_type> inline Array3D<data_type>& Array3D<data_type>::operator= (const Array3D<data_type>& a)
{
	this->resize(a.getRows(), a.getCols(), a.getSub());

      for (array_index i=0; i<a.getRows(); i++)
      {
            for (array_index j=0; j<a.getCols(); j++)
            {
                  for (array_index m=0; m<a.getSub(); m++)
                        (*this)(i,j,m) = a(i,j,m);
            }
      }

	return *this;
}


/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                    PUBLIC TEMPLATE INSTANTIATION                       ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */

template class Array1D<int>;
template class Array1D<float>;
template class Array1D<double>;
template class Array1D<string>;


template class Array2D<int>;
template class Array2D<float>;
template class Array2D<double>;

template class Array3D<int>;
template class Array3D<float>;
template class Array3D<double>;



/*
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 *  ------------------------------------------------------------------------------
 *  ###                    PRIVATE FUNCTION DEFINITIONS                        ###
 *  ------------------------------------------------------------------------------
 *  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 */




