
      program CheckDers
c
c *********************************************************************
c ******    a simple program to numerically check the rate       ******
c ******    derivatives w.r.t. temperature                       ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May       2, 2009                        ******
c ****** last updated:  March    10, 2010                        ******
c *********************************************************************
c
      Use finiteRateSurface
      Implicit none

      integer       :: ierr, i
      character*40  :: fname
      type(finiteRateSurfaceModelType) :: smodel
      real*8  :: tg1, tg2, tgmid, tw1, tw2, twmid,
     &           wm, stot, sitep, A, B, E, v, Ed, Ea
      real*8 :: frate(5,0:5), dGfrate(5,0:5), dWfrate(5,0:5),
     &          eq(5,0:4,2), dGeq(5,0:4,2), dWeq(5,0:4,2)

c     ----- SET RANDOM DUMMY VALUES -----
      wm    = 16.0d0
      stot  = 7.5d-6
      sitep = 1.0
	A     = 0.01
      B     = 1.0
      E     = 1.0d4
      
      v     = 1.0d12
      Ed    = 3.0d4
      Ea    = 4.0d4

c     ----- SET TEMPERATURES -----
      tg1 = 1000.0d0
      tg2 = 1010.0d0
      tw1 = 500.0d0
      tw2 = 501.0d0
      tgmid = 0.5*(tg1 + tg2)
      twmid = 0.5*(tw1 + tw2)

c     ----- FORWARD RATE CONTSTANT -----
      do i = 0,5
         call srate(twmid, tg1, wm, stot, sitep, i, A, B, E, 
     &              frate(1,i), dGfrate(1,i), dWfrate(1,i))
         call srate(twmid, tg2, wm, stot, sitep,  i, A, B, E, 
     &              frate(2,i), dGfrate(2,i), dWfrate(2,i))
         call srate(tw1, tgmid, wm, stot, sitep,  i, A, B, E, 
     &              frate(3,i), dGfrate(3,i), dWfrate(3,i))
         call srate(tw2, tgmid, wm, stot, sitep,  i, A, B, E, 
     &              frate(4,i), dGfrate(4,i), dWfrate(4,i))
         call srate(twmid, tgmid, wm, stot, sitep,  i, A, B, E, 
     &              frate(5,i), dGfrate(5,i), dWfrate(5,i))
      enddo

      do i = 0,5
         write(6,100) i
         write(6,101) 
         write(6,102) twmid,tgmid,frate(5,i),dGfrate(5,i),dWfrate(5,i)
         write(6,103) (frate(2,i)-frate(1,i))/(tg2-tg1),
     &                (frate(4,i)-frate(3,i))/(tw2-tw1)
         write(6,*)
         write(6,*)
      enddo
      
      write(6,300)

c     NOTE: Equilibrium Constant Type 1 is a function of forward rate.
c     therefore, forward rate derivates are needed.  The forward rate is 
c     passed at the correct temperature, but the derivatives are given
c     as computed at the average midpoints to be exact.
      
c     ----- EQUILIBRIUM CONSTANTS, TYPE 1 -----
      do i = 0,3

         call adeq(twmid, tg1, wm, 1, i, A, B, v, Ed, Ea,
     &             frate(1,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(1,i,1), dGeq(1,i,1), dWeq(1,i,1) )
         call adeq(twmid, tg2, wm, 1, i, A, B, v, Ed, Ea,
     &             frate(2,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(2,i,1), dGeq(2,i,1), dWeq(2,i,1) )
         call adeq(tw1, tgmid, wm, 1, i, A, B, v, Ed, Ea,
     &             frate(3,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(3,i,1), dGeq(3,i,1), dWeq(3,i,1) )
         call adeq(tw2, tgmid, wm, 1, i, A, B, v, Ed, Ea,
     &             frate(4,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(4,i,1), dGeq(4,i,1), dWeq(4,i,1) )
         call adeq(twmid, tgmid, wm, 1, i, A, B, v, Ed, Ea,
     &             frate(5,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(5,i,1), dGeq(5,i,1), dWeq(5,i,1) )
      enddo

      do i = 0,3
         write(6,200) 1, i
         write(6,201) 
         write(6,202) twmid,tgmid,eq(5,i,1),dGeq(5,i,1),dWeq(5,i,1)
         write(6,203) (eq(2,i,1)-eq(1,i,1))/(tg2-tg1),
     &                (eq(4,i,1)-eq(3,i,1))/(tw2-tw1)
         write(6,*)
         write(6,*)
      enddo

      write(6,300)


c     ----- EQUILIBRIUM CONSTANTS, TYPE 2 -----
      do i = 0,4

         call adeq(twmid, tg1, wm, 2, i, A, B, v, Ed, Ea,
     &             frate(1,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(1,i,2), dGeq(1,i,2), dWeq(1,i,2) )
         call adeq(twmid, tg2, wm, 2, i, A, B, v, Ed, Ea,
     &             frate(2,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(2,i,2), dGeq(2,i,2), dWeq(2,i,2) )
         call adeq(tw1, tgmid, wm, 2, i, A, B, v, Ed, Ea,
     &             frate(3,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(3,i,2), dGeq(3,i,2), dWeq(3,i,2) )
         call adeq(tw2, tgmid, wm, 2, i, A, B, v, Ed, Ea,
     &             frate(4,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(4,i,2), dGeq(4,i,2), dWeq(4,i,2) )
         call adeq(twmid, tgmid, wm, 2, i, A, B, v, Ed, Ea,
     &             frate(5,0), dGfrate(5,0), dWfrate(5,0),
     &             eq(5,i,2), dGeq(5,i,2), dWeq(5,i,2) )
      enddo

      do i = 0,4
         write(6,200) 2, i
         write(6,201) 
         write(6,202) twmid,tgmid,eq(5,i,2),dGeq(5,i,2),dWeq(5,i,2)
         write(6,203) (eq(2,i,2)-eq(1,i,2))/(tg2-tg1),
     &                (eq(4,i,2)-eq(3,i,2))/(tw2-tw1)
         write(6,*)
         write(6,*)
      enddo



 
 100  format('Forward Rate Data for Type ',i1,':')
 101  format('      Tw       Tg        frate       dF/dTg',
     &       '       dF/dTw')
 102  format(f8.1,x,f8.1,x,es12.5,x,es12.5,x,es12.5)
 103  format('   computed dF/dT values = ',4x,es12.5,x,es12.5)

 200  format('Equilibrium Constant Data for Form ',i1,' Type ',i1,':')
 201  format('      Tw       Tg           eq      dEQ/dTg',
     &       '      dEQ/dTw')
 202  format(f8.1,x,f8.1,x,es12.5,x,es12.5,x,es12.5)
 203  format('   computed dEQ/dT values = ',3x,es12.5,x,es12.5)

 300  format('=======================================================')

      end program CheckDers
