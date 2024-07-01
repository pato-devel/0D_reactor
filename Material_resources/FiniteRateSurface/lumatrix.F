
c --------------------------------------------------------------------- 
c --------------------Subroutine rowechelon()--------------------------
c ---------------------------------------------------------------------
      Subroutine rowechelon(a, irow, icol, ibasis)
c
c *********************************************************************
c ******    put matrix A in row-echelon form to find basis       ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      August   17, 2010                        ******
c ****** last updated:  August   17, 2010                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     irow   = number of rows
c     icol   = number of columns
c     a      = matrix (is destroyed)
c     ibasis = list of basis vector for each row
      Integer                     :: irow, icol
      real*8,dimension(irow,icol) :: a
      Integer,dimension(irow)     :: ibasis
c 
c *** declaration of internally used variables
c     ir     = row counter
c     indexv = ordering array
c     ans    = solution vector
      Integer :: ir, ic, ibase_col, ibase_row, iswap, iflag
      real*8  :: tmp, tgt
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 	
      ibasis = 0

      ibase_col = 0
      do ibase_row = 1,irow

         ibase_col = ibase_col + 1
c        --- search to find a base element for this row ---
         iflag = 0
         do while ((iflag .eq. 0).and.(ibase_col.le.icol))
c           --- find largest remaining row in current column ---
            tgt = 0.0d0
            iswap = ibase_row
            do ir = ibase_row, irow
               if ( abs(a(ir,ibase_col)) .gt. tgt) then
                  iswap = ir
                  tgt = abs(a(ir,ibase_col))
               endif
            enddo
c           --- determine if current column will work ---
            if (tgt .ne. 0.0d0) then
c              --- swap rows if necessary ---
               if (iswap .ne. ibase_row) then
                  do ic = 1,icol
                     tmp = a(iswap,ic)
                     a(iswap,ic)     = a(ibase_row,ic)
                     a(ibase_row,ic) = tmp
                  enddo
               endif
c              --- normalize ---
               tmp = 1.0d0/a(ibase_row,ibase_col)
               do ic = 1,icol
                  a(ibase_row,ic) = tmp*a(ibase_row,ic)
               enddo
c              --- remove remaining elements in lower rows ---
               do ir = ibase_row+1, irow
                  do ic = 1,icol
                     a(ir,ic) = a(ir,ic) - 
     &                          a(ir,ibase_col)*a(ibase_row,ic)
                  enddo
               enddo
c              --- store which column vector is basis for this row ---
               ibasis(ibase_row) = ibase_col
c              --- update flag ---
               iflag = 1
            else
               ibase_col = ibase_col + 1
            endif
         enddo    !do-while loop
         
      enddo   !ibase_row loop
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine rowechelon
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c -----------------------Subroutine lusolve()--------------------------
c ---------------------------------------------------------------------
      Subroutine lusolve(alhs, brhs, isize)
c
c *********************************************************************
c ******    solve the system Ax = b                              ******
c ******    store the result in b; A is destroyed                ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      18, 2009                        ******
c ****** last updated:  May      18, 2009                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     isize = dimension of square matrix
c     alhs  = left hand side matrix
c     brhs  = right hand side matrix (overridden)
      Integer                       :: isize
      real*8,dimension(isize,isize) :: alhs
      real*8,dimension(isize)       :: brhs
c 
c *** declaration of internally used variables
c     ir     = row counter
c     indexv = ordering array
c     ans    = solution vector
      Integer                       :: ir
      real*8,dimension(isize)       :: ans
      integer,dimension(isize)      :: indexv
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c     ----- compute LU matrix -----
      call ludecomp(alhs, indexv, isize)

c     ----- reorder RHS matrix -----      
      call lupivotrhs(brhs, indexv, isize)

c     ----- solve for solution vector -----      
      call lusubst(alhs, brhs, ans, isize)

c     -----substitute into RHS matrix -----
      do ir = 1,isize
         brhs(ir) = ans(ir)
      enddo
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine lusolve
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ---------------------Subroutine luinverse()--------------------------
c ---------------------------------------------------------------------
      Subroutine luinverse(amat, isize)
c
c *********************************************************************
c ******    invert a matrix using LU decomposition               ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      18, 2009                        ******
c ****** last updated:  May      18, 2009                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     size  = dimension of square matrix
c     alu   = matrix to invert (overridden)
      Integer                       :: isize
      real*8,dimension(isize,isize) :: amat
c 
c *** declaration of internally used variables
c     ir  = row counter
c     tmp = ordering array
      Integer                       :: ir, ic, k, iv
      real*8,dimension(isize)       :: rhs, x
      integer,dimension(isize)      :: indexv
      real*8,dimension(isize,isize) :: atmp
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      atmp = 0.0d0

c     ----- decompose the matrix A once -----
      call ludecomp(amat, indexv, isize)

c     ----- find columns of the inverse by successive substitution -----
      do ic = 1,isize

         rhs     = 0.0d0
         rhs(ic) = 1.0d0
         call lusubst(amat, rhs, x, isize)
         
         iv = indexv(ic)
         do ir = 1,isize
            atmp(ir,iv) = x(ir)
         enddo
         
      enddo

c     ----- overwrite the passed memory block -----      
      do ic = 1,isize
         do ir = 1,isize
            amat(ir,ic) = atmp(ir,ic)
         enddo
      enddo
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine luinverse
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ---------------------Subroutine ludecomp()---------------------------
c ---------------------------------------------------------------------
      Subroutine ludecomp(amat, indexv, isize)
c
c *********************************************************************
c ******    decompose matrix A into lower and upper symmetric    ******
c ******    components with pivoting reordering                  ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      16, 2009                        ******
c ****** last updated:  May      16, 2009                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     size     = dimension of square matrix
c     amat     = square matrix to invert
c     indexv   = index reordering vector
      Integer                       :: isize
      real*8,dimension(isize,isize) :: amat
      integer,dimension(isize)      :: indexv
c 
c **** declaration of internally used variables
c      i,n         = loop counters
c      nrange      = loop counters for bulk and surface indices
c      imatch      = index of a matched species
c      ier         = internal error for file problems
c      coefftmp    = array of coefficients read from file
c      snm         = name of species read from file
c      ifound      = array to track whether fits were 
c                    found for each species
      Integer                 :: ir, ic, k, ipindex, itmp
      real*8                  :: pivot, tmp
      real*8,dimension(isize) :: scale
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 


c     ----- initialize index order -----      
      do k = 1,isize
         indexv(k) = k
      enddo


c     ----- find largest element in each row -----
      do ir = 1,isize
         scale(ir) = 0.0d0
         do ic = 1,isize
            tmp = abs( amat(ir,ic) )
            if ( tmp .gt. scale(ir) ) scale(ir) = tmp
         enddo
         scale(ir) = 1.0d0 / scale(ir)
      enddo


c     ----- main loop -----
      do ic = 1,isize

c        ----- operate on indices up to but not including the diagonal -----
         do ir = 1,ic-1
            do k = 1,ir-1
               amat(ir,ic) = amat(ir,ic) - amat(ir,k) * amat(k,ic)
            enddo
         enddo
         
         
         pivot = 0.0d0
         do ir = ic,isize
            do k = 1,ic-1
               amat(ir,ic) = amat(ir,ic) - amat(ir,k) * amat(k,ic)
            enddo
            tmp = scale(ir) * abs( amat(ir,ic) )
            if (tmp > pivot) then
               pivot   = tmp
               ipindex = ir
            endif
         enddo
         
         
         if (ipindex .ne. ic) then
            do k = 1,isize
               tmp = amat(ipindex,k)
               amat(ipindex,k) = amat(ic,k)
               amat(ic,k) = tmp
            enddo
            itmp = indexv(ic)
            indexv(ic) = indexv(ipindex)
            indexv(ipindex) = itmp
            scale(ipindex) = scale(ic)
         endif
         
         
         do ir = ic+1,isize
            amat(ir,ic) = amat(ir,ic) / amat(ic,ic)
         enddo
         
      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine ludecomp
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c -------------------Subroutine lupivotrhs()---------------------------
c ---------------------------------------------------------------------
      Subroutine lupivotrhs(rhs, indexv, isize)
c
c *********************************************************************
c ******    pivot a rhs vector according to index order          ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      16, 2009                        ******
c ****** last updated:  May      16, 2009                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     size     = dimension of square matrix
c     amat     = square matrix to invert
c     indexv   = index reordering vector
      Integer                   :: isize
      real*8,dimension(isize)   :: rhs
      integer,dimension(isize)  :: indexv
c 
c *** declaration of internally used variables
c     ir  = row counter
c     tmp = ordering array
      Integer                 :: ir
      real*8,dimension(isize) :: tmp
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

      do ir = 1,isize
         tmp(ir) = rhs(indexv(ir))
      enddo
      
      do ir = 1,isize
         rhs(ir) = tmp(ir)
      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine lupivotrhs
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ---------------------Subroutine lusubst()----------------------------
c ---------------------------------------------------------------------
      Subroutine lusubst(alu, rhs, dx, isize)
c
c *********************************************************************
c ******    solve a system that has undergone LU decomposition   ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      16, 2009                        ******
c ****** last updated:  May      16, 2009                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     size  = dimension of square matrix
c     alu   = LU stored matrix
c     rhs   = rhs vector (already pivoted)
c     dx    = solution vector
      Integer                       :: isize
      real*8,dimension(isize)       :: rhs, dx
      real*8,dimension(isize,isize) :: alu
c 
c *** declaration of internally used variables
c     ir  = row counter
c     tmp = ordering array
      Integer                 :: ir, k
      real*8,dimension(isize) :: tmp
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c     ----- forward substitution step -----
      do ir = 1,isize
         tmp(ir) = rhs(ir)
         do k=1,ir-1
            tmp(ir) = tmp(ir) - alu(ir,k)*tmp(k)
         enddo
      enddo

c     ----- backward substitution step -----
      do ir = isize,1,-1
         dx(ir) = tmp(ir)
         do k=ir+1,isize
            dx(ir) = dx(ir) - alu(ir,k)*dx(k)
         enddo
         dx(ir) = dx(ir) / alu(ir,ir)
      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine lusubst
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
