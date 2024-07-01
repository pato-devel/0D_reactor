
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      program SRmodel
c
c     J. Marschall
c     2/2009
c
c     This program reads surface reaction model inputs,
c     calculates and stores reaction rates and equilibrium
c     constants, and assembles equations for surface coverage
c  
c     Not pretty FORTRAN
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Use finiteRateSurface
      Implicit none
c
c *** parameters
      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
      
      integer :: ierr, i, k, ngas, iphase, ireact, ispec
      real*8  :: tw, tg, p, tmp, prate, pf, pb, vave
      character*40  :: gname
      character*128 :: fname
	character*200 :: bname
      real*8,dimension(100)  :: frates,dfg,dfw,eqa,eqc,brates,gloss,
     &                          rflux
      real*8,dimension(8)    :: gibbs,dgg,dgw
      
      character*6,dimension(5) :: gspecies
      real*8,dimension(10,1,5) :: gcof
      real*8,dimension(2)      :: gtmpr
      real*8,dimension(5)      :: gmolefr
      
      type(finiteRateSurfaceModelType) :: frmr
      type(finiteRateSurfaceCell)      :: cell

c     ----- set gas phase species names -----
      gspecies(1) = 'N2'
      gspecies(2) = 'O2'
      gspecies(3) = 'NO'
      gspecies(4) = 'N'
      gspecies(5) = 'O'
      tw = 2000.
      tg = 2000.
      p = 10000.
      
      gmolefr(1) = 0.690
      gmolefr(2) = 0.150
      gmolefr(3) = 0.010
      gmolefr(4) = 0.050
      gmolefr(5) = 0.100
      

c     ----- read in surface model file -----
      fname = 'surf/air_silica.inp'
      call read_finite_rate_surf_file(fname, frmr, ierr)

      
c     ----- set fake Gibbs energies to match Joe's data -----
c           need to pass non-dimensional data to Lewis 
c           routine
      gcof = 0.0d0
      do i = 1,5
         gcof(10,1,i) = -100000.0d0/(R*tg)
      enddo
      gtmpr(1) = 0.0d0
      gtmpr(2) = 1.0d10

c     -----      
      bname = 'lewis.bulk'
      call read_lewis_bulk_file(bname, frmr, gcof, gtmpr, 1, ierr)


c---- initialize the cell object
      call finite_rate_surf_initialize(cell, frmr)
      cell%gas_tmpr  = tg
      cell%wall_tmpr = tw
      cell%wall_pres = p
      do k=1,5
        cell%concentration(k) = p*gmolefr(k)/(R*tg)
      enddo


c---- work through all forward reaction rates one at a time
      do ireact = 1,frmr%isrt

c        ---- identify gas phase participant ----      
         k = frmr%RR(ireact,1)
         if (frmr%istype(ireact) .eq. 4) k = frmr%RP(ireact,1)

c        ---- compute forward rate and derivatives ----
         call srate(tw, tg, frmr%wmass(k),
     &              frmr%sdensity(ireact), frmr%sitep(ireact),
     &              frmr%istype(ireact),
     &              frmr%Cf(ireact), frmr%eta(ireact),
     &              frmr%Eact(ireact), frates(ireact),
     &              dfg(ireact), dfw(ireact))
   
      enddo


c---- obtain Gibbs energies for all species either by Lewis data or backward rates
      call finite_rate_surf_gibbs_energy(frmr, gibbs, dgg, dgw,
     &                                   tg, tw, frates, dfg, dfw)


c---- calculate activity based equilibrium constants (needed only for output)
      do ireact = 1,frmr%isrt
         eqa(ireact) = 0.0d0
         do ispec = 1,frmr%nspecies
            tmp = frmr%v(ireact,ispec)/(R*tw)
            eqa(ireact) = eqa(ireact) - tmp*gibbs(ispec)
         enddo
         eqa(ireact) = exp(eqa(ireact))
      enddo


c---- calculate concentration based equilibrium constants
      call finite_rate_surf_eqcon(frmr, gibbs, dgg, dgw, 
     &                            cell%wall_tmpr, eqc, dfg, dfw)


c---- calculate backward rates
      do ireact = 1,frmr%isrt
         brates(ireact) = frates(ireact)/eqc(ireact)
      enddo


c---- calculate flux for each reaction
      do ireact = 1,frmr%isrt

c        ----- sum forward and backward concentration product term -----
         pf = 1.0d0
         pb = 1.0d0
         do ispec = 1,frmr%nspecies
            pf = pf*cell%concentration(ispec)**frmr%vr(ireact,ispec)
            pb = pb*cell%concentration(ispec)**frmr%vp(ireact,ispec)
         enddo

c        ---- account for rates turned on and compute reaction flux ----
         pf = pf * frmr%isrfon(ireact)
         pb = pb * frmr%isrbon(ireact)
         rflux(ireact) = frates(ireact)*pf - brates(ireact)*pb

c        ---- compute this contribution to species production ----         
         do ispec = 1,frmr%nspecies
            prate = frmr%v(ireact,ispec) * rflux(ireact)

c           ---- explicit contribution ----           
            cell%srhs(ispec) = cell%srhs(ispec) + prate

         enddo
         
      enddo

c---- calculate loss efficiency for gas phase
      do k=1,frmr%ngps
        vave=sqrt(8.0*R*Tg/(PI*frmr%wmass(k)) )
        gloss(k)=-cell%srhs(k)/(cell%concentration(k)*vave/4.0)
      enddo      



c-------------------------------------------------
c-----Read in all phase and species information
c-------------------------------------------------

c---- Write info
      write(6,*) " "
      write(6,*) "Species  Species#    Molar Mass    Phase     Phase# "
      do k=1,frmr%nspecies 
        write(6,111) frmr%species(k),k,frmr%wmass(k),
     &               frmr%phases(frmr%kphase(k)),
     &               frmr%kphase(k)
      enddo

c-----------------------------------------------
c-----Read in reaction information--------------
c-----------------------------------------------

c---- Write info.
      write(6,*) " "
      write(6,*) "Reactant and product identification for each reaction"
      do i=1,frmr%isrt
        write(6,*) frmr%RR(i,1),frmr%RR(i,2),frmr%RR(i,3),
     &             frmr%RP(i,1),frmr%RP(i,2),frmr%RP(i,3)
      enddo
      write(6,*) " "
      write(6,*) "Reactant and product stoichiometry for each reaction"
      do i=1,frmr%isrt
        write(6,*) frmr%RVR(i,1),frmr%RVR(i,2),frmr%RVR(i,3),
     &             frmr%RVP(i,1),frmr%RVP(i,2),frmr%RVP(i,3)
      enddo
      write(6,*) " "
      write(6,*) "Reactant stoichiometry matrix"
      do i=1,frmr%isrt
        write(6,122) (frmr%vr(i,k),k=1,frmr%nspecies )
      enddo
      write(6,*) " "
      write(6,*) "Product stoichiometry matrix"
      do i=1,frmr%isrt
        write(6,122) (frmr%vp(i,k),k=1,frmr%nspecies )
      enddo
      write(6,*) " "
      write(6,*) "Product minus Reactant stoichiometry matrix"
      do i=1,frmr%isrt
        write(6,122) (frmr%v(i,k),k=1,frmr%nspecies )
      enddo
      write(6,*) " "

c-----------------------------------------------------
c---- Read in rate constant information --------------
c-----------------------------------------------------

c---- Write info.
      write(6,*) " "
      write(6,*) "Reaction Type  Cf    eta  Ea   isrfon   isrbon"
      do i=1,frmr%isrt
        write(6,*) i,frmr%istype(i),frmr%Cf(i),frmr%eta(i),
     &               frmr%Ea(i),frmr%isrfon(i),frmr%isrbon(i)
      enddo

c--------------------------------------------------------------
c---- For surface phases without thermochemical data ----------
c---- Read in desorption reactions or equilibrium constants ---
c--------------------------------------------------------------

c---- Write info.
      write(6,*) " "
      write(6,*) "Reaction Type  Form  Cf  eta  vdes  Edes"
      do i=1,frmr%isrt
        write(6,*) i,frmr%iadtype(i),frmr%iadform(i),
     &             frmr%Cfad(i),frmr%etaad(i),
     &             frmr%vad(i),frmr%Edes(i)
      enddo
      write(6,*) " "


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c---- The following would be repeated for each set of surface temperatures
c---- and species concentrations
c-----------------------------------------------------------------------
c------------------------------------------------------------------------
c---- Write info.
      write(6,*) " "
      write(6,*) "Species   Generalized Concentration "
      do k=1,frmr%nspecies
        write(6,*) k,frmr%species(k),cell%concentration(k)
      enddo
      write(6,*) " "
        
               
c--------------------------------------------
c---- Evaluate forward rate coefficients ----
c--------------------------------------------
c---- Write info.
      write(6,*) " "
      write(6,*) "Reaction  Type  Rate "
      do i=1,frmr%isrt
        write(6,*) i,frmr%istype(i),frates(i)
      enddo
      write(6,*) " "


c--------------------------------------------
c---- Collect Gibbs Formation Energy  Data --
c--------------------------------------------

c---- write info.
      write(6,*) " "
      write(6,*) "Species# Phase# ThermoData  AdsReac   GibbsFE "
      do k=1,frmr%nspecies 
        write(6,*) k,frmr%kphase(k),frmr%kthermo(k),
     &             frmr%kadsr(k),gibbs(k)
      enddo
      write(6,*) " "


c-----------------------------------------------------------
c---- Calculate concetration-based equilibrium constants ---
c-----------------------------------------------------------

c---- Write info
      write(6,*) "Reaction  Kqa  Kqc"
      do i=1,frmr%isrt
        write(6,*) i,eqa(i),eqc(i)
      enddo
      write(6,*) " "


c---------------------------------------------------------
c---- Calculate reaction flux for each reaction ----------
c---------------------------------------------------------

c---- Write info
      write(6,*) "Reaction  frate  fon  brate  bon  rflux"
      do i=1,frmr%isrt
        write(6,*) i,frates(i),frmr%isrfon(i),
     &               brates(i),frmr%isrbon(i),rflux(i)
      enddo
      write(6,*) " "


c--------------------------------------------------------------
c---- Calculate the loss efficiency for each gas species ------
c--------------------------------------------------------------

c---- Write info
      write(6,*) "Species   production rate   loss efficiency"
      do k=1,frmr%ngps
         write(6,*) frmr%species(k),cell%srhs(k),gloss(k)
      enddo
      do k=frmr%ngps+1,frmr%nspecies
         write(6,*) frmr%species(k),cell%srhs(k)
      enddo


 111   format(1x,a10,i5,5x,f10.8,6x,a10,i5)
 122   format(1x,21i3)

      end program
