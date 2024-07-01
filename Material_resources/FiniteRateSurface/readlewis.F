c --------------------------------------------------------------------- 
c ----------------Subroutine read_lewis_gas_file()---------------------
c ---------------------------------------------------------------------
      Subroutine read_lewis_gas_file(fname, nspecies, snames, 
     &                               cof, ierr)
c
c *********************************************************************
c ******    DPLR provides the gas phase Lewis data, but we       ******
c ******    need this routine to run the offline test code(s)    ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May       8, 2009                        ******
c ****** last updated:  May       2, 2011                        ******
c *********************************************************************
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     fname    = name of file to open
c     snames   = array of species names to look for
c     nspecies = number of species in search
c     ierr     = error flag for return
c     cof      = gas phase Gibbs coefficient array to fill-in
      Character*40                           :: fname
      Character*6,dimension(nspecies)        :: snames
      type(finiteRateSurfaceModelType)       :: smodel
      Integer                                :: ierr, nspecies
      real*8,dimension(10,5,nspecies)        :: cof
c 
c **** declaration of internally used variables
c      i,j,n       = loop counters
c      m           = range index
c      nrange      = loop counters for bulk and surface indices
c      imatch      = index of a matched species
c      ier         = internal error for file problems
c      coefftmp    = array of coefficients read from file
c      snm         = name of species read from file
c      h,s,cp      = Lewis functions
c      t(*)        = temperature functions for Lewis evaluations
c      ifound      = array to track whether fits were 
c                    found for each species
      Integer                       :: nrange,imatch,i,j,m,n,ier
      real*8                        :: coefftmp(10),
     &                                 t1,t2,t3,t4,ti,ts,tl,h,s,cp
      Character*8                   :: snm
      Integer,dimension(0:nspecies) :: ifound
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

      coefftmp = 0.0d0
      cof      = 0.0d0
      ifound   = 0
      ierr = 0


      open(10, file=fname, status='old',iostat=ier)
      if (ier .ne. 0) then
        ierr = -1
        return
      endif
      rewind(10)
      read(10,*) 


c---------------------------------------------------------
c-----pass through file looking for matching species names
c---------------------------------------------------------

      do j = 1,1000

         read(10,*) 
         read(10,*) snm,nrange

         if (trim(snm) .eq. 'END') goto 121

c     ----- look for a gas species that matches this name -----
         imatch = 0
         do n = 1,nspecies
            if ( trim(snm) .eq. trim(snames(n)) ) imatch = n
         enddo
         ifound(imatch) = 1

         do i = 1,nrange
            read(10,*)  coefftmp(1:5)
            read(10,*)  coefftmp(6:10)
            if (imatch .gt. 0) then
               do n = 1,10
                  cof(n,i+1,imatch) = coefftmp(n)
               enddo
            endif
         enddo

c     ----- extend ranges on Lewis fits assuming constant Cp -----
         if (imatch .gt. 0) then

c           ----- first range -----
            m = 2
            t1 = 200.0d0
            t2 = t1*t1
            t3 = t2*t1
            t4 = t3*t1
            ti = 1.0d0/t1
            ts = ti*ti
            tl = log(t1)
            h  = cof(3,m,imatch)
     &         + cof(2,m,imatch)*ti*tl
     &         - cof(1,m,imatch)*ts
     &         + 0.5*cof(4,m,imatch)*t1
     &         + 1.0d0/3.0d0*cof(5,m,imatch)*t2
     &         + 0.25*cof(6,m,imatch)*t3
     &         + 0.2*cof(7,m,imatch)*t4
     &         + cof(9,m,imatch)*ti
            s  = cof(3,m,imatch)*tl
     &         - cof(2,m,imatch)*ti
     &         - 0.5*cof(1,m,imatch)*ts
     &         + cof(4,m,imatch)*t1
     &         + 0.5*cof(5,m,imatch)*t2
     &         + 1.0d0/3.0d0*cof(6,m,imatch)*t3
     &         + 0.25*cof(7,m,imatch)*t4
     &         + cof(10,m,imatch)
            cp = cof(1,m,imatch)*ts
     &         + cof(2,m,imatch)*ti
     &         + cof(3,m,imatch)
     &         + cof(4,m,imatch)*t1
     &         + cof(5,m,imatch)*t2
     &         + cof(6,m,imatch)*t3
     &         + cof(7,m,imatch)*t4

            cof(3,m-1,imatch)  = cp
            cof(9,m-1,imatch)  = (h-cp)*t1
            cof(10,m-1,imatch) = s - cp*tl

c           ----- end range -----
            m = nrange + 1
            if (nrange .eq. 1) t1 = 1000.0d0
            if (nrange .eq. 2) t1 = 6000.0d0
            if (nrange .eq. 3) t1 = 20000.0d0
            t2 = t1*t1
            t3 = t2*t1
            t4 = t3*t1
            ti = 1.0d0/t1
            ts = ti*ti
            tl = log(t1)
            h  = cof(3,m,imatch)
     &         + cof(2,m,imatch)*ti*tl
     &         - cof(1,m,imatch)*ts
     &         + 0.5*cof(4,m,imatch)*t1
     &         + 1.0d0/3.0d0*cof(5,m,imatch)*t2
     &         + 0.25*cof(6,m,imatch)*t3
     &         + 0.2*cof(7,m,imatch)*t4
     &         + cof(9,m,imatch)*ti
            s  = cof(3,m,imatch)*tl
     &         - cof(2,m,imatch)*ti
     &         - 0.5*cof(1,m,imatch)*ts
     &         + cof(4,m,imatch)*t1
     &         + 0.5*cof(5,m,imatch)*t2
     &         + 1.0d0/3.0d0*cof(6,m,imatch)*t3
     &         + 0.25*cof(7,m,imatch)*t4
     &         + cof(10,m,imatch)
            cp = cof(1,m,imatch)*ts
     &         + cof(2,m,imatch)*ti
     &         + cof(3,m,imatch)
     &         + cof(4,m,imatch)*t1
     &         + cof(5,m,imatch)*t2
     &         + cof(6,m,imatch)*t3
     &         + cof(7,m,imatch)*t4

            do i=(m+1),5
               cof(3,i,imatch)  = cp
               cof(9,i,imatch)  = (h-cp)*t1
               cof(10,i,imatch) = s - cp*tl
            enddo

         endif

      enddo

c     ----- close file -----
 121  close (10)

c     ----- make sure all species found or issue error -----
      do n = 1,nspecies
         if ( ifound(n) .eq. 0) ierr = n
      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine read_lewis_gas_file
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

