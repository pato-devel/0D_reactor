
      program CheckSystem
c
c *********************************************************************
c ******    a program to numerically check all system derivative ******
c ******    components against the analytical values             ******
c ******                                                         ******
c ******    This is essentially just a diagnostic version        ******
c ******    of the finite_rate_compute_system subroutine         ******
c ******    located in FiniteRateSurface.F                       ******
c ******                                                         ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      May      12, 2009                        ******
c ****** last updated:  April    15, 2011                        ******
c *********************************************************************
c
      Use finiteRateSurface
      Implicit none

      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
      
      real*8 ::DTMPR, DCONC
      parameter(DTMPR=2.0d0,DCONC=0.01)

      integer :: ierr, i, k, ngas, it, ireact, iphase, ispec, nsteps,
     &           iconv, iout, itorder, idiagnostic, ieqtype,
     &           igastype, iscan
      real*8  :: tw(5), tg(5), p, dt, smoltot
      character*40  :: gname, pname
	character*128 :: fname
	character*200 :: bname
      character*6,allocatable,dimension(:) :: snames
      real*8,allocatable,dimension(:)      :: smolefr, coftmpr
      real*8,allocatable,dimension(:,:,:)  :: cof
      
      real*8,allocatable,dimension(:,:) :: frates, brates, 
     &                                     eqc, gibbs, rhs, rnumjac
      real*8,allocatable,dimension(:)   :: ftmp, eqtmp, gtmp,
     &                                     dfg, dfw, dbg, dbw,
     &                                     dcg, dcw, dgg, dgw
      
      type(finiteRateSurfaceModelType) :: frmr
      type(finiteRateSurfaceCell)      :: cell


c-------------------------------------------------
c-----Read in input file
c-------------------------------------------------

      open(13, file='input.inp', status='old')
      
      read (13,*)
      read (13,*) fname
      
      read (13,*)
      read (13,*)
      read (13,*) iconv, iout, pname
      
      read (13,*)
      read (13,*)
      read (13,*) tw(5), p, igastype, iscan
      tg(5) = tw(5)

      read (13,*)
      read (13,*)
      read (13,*) dt, nsteps, itorder, idiagnostic, ieqtype
      
      read (13,*)
      read (13,*)
      read (13,*) ngas
      
      allocate(snames(ngas),smolefr(ngas))
      
      read (13,*)
      read (13,*)
      do i = 1,ngas
         read (13,*) snames(i), smolefr(i)
      enddo

      close(13)
      
c-------------------------------------------------
c-----Read in input files and collect data
c-------------------------------------------------

c     ----- read in data from surface file -----
      call read_finite_rate_surf_file(fname, frmr, ierr)

c     ----- read gas phase lewis/gibbs data -----      
      allocate(cof(10,5,ngas), coftmpr(6))
      gname = 'lewis.thermo'
      call read_lewis_gas_file(gname, ngas, snames, 
     &                         cof, ierr)
      coftmpr(1) = 0.0d0
      coftmpr(2) = 200.0d0
      coftmpr(3) = 1000.0d0
      coftmpr(4) = 6000.0d0
      coftmpr(5) = 20000.0d0
      coftmpr(6) = 1.0d+20

c     ----- read bulk phase lewis/gibbs data -----
      bname = 'lewis.bulk'
      call read_lewis_bulk_file(bname, frmr, cof, coftmpr, 3, ierr)


c-------------------------------------------------
c-----Initialize single surface cell
c-------------------------------------------------

c     ----- initialize and allocate arrays -----
      call finite_rate_surf_initialize(cell, frmr)

c     ----- set values by hand -----
      cell%gas_tmpr  = tg(5)
      cell%wall_tmpr = tw(5)
      cell%wall_pres = p

c     ----- normalize mole fraction values -----
      smoltot = 0.0d0
      do i = 1,ngas
         if (smolefr(i) .lt. 1.0d-15) smolefr(i) = 1.0d-15
         smoltot = smoltot + smolefr(i)
      enddo
      do i = 1,ngas
         smolefr(i) = smolefr(i)/smoltot
      enddo

c     ----- override gas phase concentrations -----
      do k=1,ngas
        cell%concentration(k) = cell%wall_pres *
     &                                    smolefr(k)/(R*cell%gas_tmpr)
      enddo


c----------------------------------------------------------
c-----Declare individual arrays & temperature perturbations
c----------------------------------------------------------

      allocate(frates(frmr%isrt,5),brates(frmr%isrt,5),
     &         eqc(frmr%isrt,5),gibbs(frmr%nspecies,5)    )
      allocate(ftmp(frmr%isrt), eqtmp(frmr%isrt), gtmp(frmr%nspecies) )
      allocate(dfg(frmr%isrt), dfw(frmr%isrt),
     &         dbg(frmr%isrt), dbw(frmr%isrt),
     &         dcg(frmr%isrt), dcw(frmr%isrt),
     &         dgg(frmr%nspecies), dgw(frmr%nspecies) )

      tg(1) = tg(5) - DTMPR
      tg(2) = tg(5) + DTMPR
      tg(3) = tg(5)
      tg(4) = tg(5)
      tw(1) = tw(5)
      tw(2) = tw(5)
      tw(3) = tw(5) - DTMPR
      tw(4) = tw(5) + DTMPR


c############################################################
c#  Run through all steps to assemble the production rates.
c#  Note that we need to compute all explicit values 5 times
c#  so that we have temperature perturbations around the 
c#  means, but we need only track the derivatives at the 
c#  mean itself.
c############################################################

c-------------------------------------------------
c-----Compute forward rates
c-------------------------------------------------

      do it = 1,5
         cell%wall_tmpr = tw(it)
         cell%gas_tmpr  = tg(it)
         do ireact = 1,frmr%isrt

c           ---- identify gas phase participant ----      
            k = frmr%RR(ireact,1)
            if (frmr%istype(ireact) .eq. 4) k = frmr%RP(ireact,1)

c           ---- compute forward rate and derivatives ----
            call srate(cell%wall_tmpr, cell%gas_tmpr, frmr%wmass(k),
     &                 frmr%sdensity(ireact), frmr%sitep(ireact),
     &                 frmr%istype(ireact),
     &                 frmr%Cf(ireact), frmr%eta(ireact),
     &                 frmr%Eact(ireact), frates(ireact, it),
     &                 dfg(ireact), dfw(ireact))
         enddo
      enddo


c-------------------------------------------------
c-----Compute Gibbs energies
c-------------------------------------------------

      do it = 1,5
         cell%wall_tmpr = tw(it)
         cell%gas_tmpr  = tg(it)
         ftmp(1:frmr%isrt) = frates(1:frmr%isrt,it)
         call finite_rate_surf_gibbs_energy(frmr, gtmp, dgg, dgw,
     &                                   cell%gas_tmpr, cell%wall_tmpr, 
     &                                   ftmp, dfg, dfw)
         gibbs(1:frmr%nspecies,it) = gtmp(1:frmr%nspecies)
      enddo



c-------------------------------------------------
c-----Compute equilibrium constants
c-------------------------------------------------

      do it = 1,5
         cell%wall_tmpr = tw(it)
         cell%gas_tmpr  = tg(it)
         ftmp(1:frmr%isrt) = frates(1:frmr%isrt,it)
         gtmp(1:frmr%nspecies) = gibbs(1:frmr%nspecies,it)
         call finite_rate_surf_eqcon(frmr, gtmp, dgg, dgw,
     &                   cell%wall_tmpr, eqtmp, dcg, dcw)
         eqc(1:frmr%isrt,it) = eqtmp(1:frmr%isrt)
      enddo


c-------------------------------------------------
c-----Compute backward rates
c-------------------------------------------------

      do it = 1,5
         do ireact = 1,frmr%isrt
            brates(ireact,it) = frates(ireact,it)/eqc(ireact,it)
         enddo
      enddo
      do ireact = 1,frmr%isrt
         dbg(ireact) = (dfg(ireact) - brates(ireact,5)*
     &                                    dcg(ireact))/eqc(ireact,5)
         dbw(ireact) = (dfw(ireact) - brates(ireact,5)*
     &                                    dcw(ireact))/eqc(ireact,5)
      enddo


c-------------------------------------------------
c-----Compute and check Jacobian for system by
c-----perturbing all concentrations a bit
c-------------------------------------------------

      allocate(rhs(frmr%nspecies,2),
     &         rnumjac(frmr%nspecies,frmr%nspecies+2) )

      rnumjac = 0.0d0
      
c      cell%concentration(7) = 0.1d-15
c      cell%concentration(8) = 0.2d-15

c     ----- compute w.r.t. each species -----
      do ispec = 1,frmr%nspecies
         
         cell%concentration(ispec) = cell%concentration(ispec) - DCONC
         call finite_rate_compute_system(cell)
         rhs(1:frmr%nspecies,1) = cell%srhs(1:frmr%nspecies)
         
         cell%concentration(ispec) =
     &                           cell%concentration(ispec) + 2.0*DCONC
         call finite_rate_compute_system(cell)
         rhs(1:frmr%nspecies,2) = cell%srhs(1:frmr%nspecies)
         
         cell%concentration(ispec) = cell%concentration(ispec) - DCONC
         
         rnumjac(1:frmr%nspecies,ispec) = 
     &      (rhs(1:frmr%nspecies,2)-rhs(1:frmr%nspecies,1))/(2.0*DCONC)
         
      enddo

c     ----- compute w.r.t. gas temperature -----
      cell%gas_tmpr = cell%gas_tmpr - DTMPR
      call finite_rate_compute_system(cell)
      rhs(1:frmr%nspecies,1) = cell%srhs(1:frmr%nspecies)
         
      cell%gas_tmpr = cell%gas_tmpr + 2.0*DTMPR
      call finite_rate_compute_system(cell)
      rhs(1:frmr%nspecies,2) = cell%srhs(1:frmr%nspecies)
         
      cell%gas_tmpr = cell%gas_tmpr - DTMPR
         
      rnumjac(1:frmr%nspecies,frmr%nspecies+1) = 
     &      (rhs(1:frmr%nspecies,2)-rhs(1:frmr%nspecies,1))/(2.0*DTMPR)

c     ----- compute w.r.t. surface temperature -----
      cell%wall_tmpr = cell%wall_tmpr - DTMPR
      call finite_rate_compute_system(cell)
      rhs(1:frmr%nspecies,1) = cell%srhs(1:frmr%nspecies)
         
      cell%wall_tmpr = cell%wall_tmpr + 2.0*DTMPR
      call finite_rate_compute_system(cell)
      rhs(1:frmr%nspecies,2) = cell%srhs(1:frmr%nspecies)

      cell%wall_tmpr = cell%wall_tmpr - DTMPR

      rnumjac(1:frmr%nspecies,frmr%nspecies+2) = 
     &      (rhs(1:frmr%nspecies,2)-rhs(1:frmr%nspecies,1))/(2.0*DTMPR)

c     ----- compute the analytic thing -----
      call finite_rate_compute_system(cell)


c-------------------------------------------------
c-----Write Output
c-------------------------------------------------

      write(6,*)
      write(6,101) 'FORWARD RATES'
      write(6,102)
      write(6,103)
      do i = 1,frmr%isrt
         write(6,201) i,frates(i,5),
     &                dfg(i),
     &                (frates(i,2)-frates(i,1))/(2.0*DTMPR),
     &                dfw(i),
     &                (frates(i,4)-frates(i,3))/(2.0*DTMPR)
      enddo


      write(6,*)
      write(6,101) 'GIBBS ENERGY'
      write(6,102)
      write(6,103)
      do i = 1,frmr%nspecies
         write(6,201) i,gibbs(i,5),
     &                dgg(i),
     &                (gibbs(i,2)-gibbs(i,1))/(2.0*DTMPR),
     &                dgw(i),
     &                (gibbs(i,4)-gibbs(i,3))/(2.0*DTMPR)
      enddo

      write(6,*)
      write(6,101) 'CONCENTRATION BASED EQUILIBRIUM CONSTANTS'
      write(6,102)
      write(6,103)
      do i = 1,frmr%isrt
         write(6,201) i,eqc(i,5),
     &                dcg(i),
     &                (eqc(i,2)-eqc(i,1))/(2.0*DTMPR),
     &                dcw(i),
     &                (eqc(i,4)-eqc(i,3))/(2.0*DTMPR)
      enddo

      write(6,*)
      write(6,101) 'BACKWARD RATES'
      write(6,102)
      write(6,103)
      do i = 1,frmr%isrt
         write(6,201) i,brates(i,5),
     &                dbg(i),
     &                (brates(i,2)-brates(i,1))/(2.0*DTMPR),
     &                dbw(i),
     &                (brates(i,4)-brates(i,3))/(2.0*DTMPR)
      enddo

      write(6,*)
      write(6,101) 'CELL RHS'
      do ispec = 1,frmr%nspecies
         write(6,301) cell%srhs(ispec)
      enddo

      write(6,*)
      write(6,101) 'ANALYTICAL JACOBIAN'
      do ispec = 1,frmr%nspecies
         write(6,301) cell%sjacobian(ispec,1:frmr%nspecies+2)
      enddo
      write(6,101) 'NUMERICAL JACOBIAN'
      do ispec = 1,frmr%nspecies
         write(6,301) rnumjac(ispec,1:frmr%nspecies+2)
      enddo

 101  format('***** ',a,' *****')
 102  format(25x,'---------Gas Tmpr--------',
     &        7x,'---------Wall Tmpr-------')
 103  format(10x,'value',12x,'analytic     numerical',10x,
     &                       'analytic     numerical')
 201  format(i3,x,es14.5,4x,2es14.5,4x,2es14.5)

 301  format(50es11.3)

      end program CheckSystem

