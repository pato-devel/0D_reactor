c --------------------------------------------------------------------- 
c -------------Subroutine finite_rate_compute_system()-----------------
c ---------------------------------------------------------------------
      Subroutine finite_rate_compute_system(cell)
c
c *********************************************************************
c ******    compute LEFT and RIGHT hand side matrices for        ******
c ******    the surface system.                                  ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ******                                                         ******
c ****** coded on:      May       9, 2009                        ******
c ****** last updated:  April    20, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c
c *** parameters
      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
c 
c *** declarations of argument types passed to the subroutine
c     cell = cell with current concentrations already set
      type(finiteRateSurfaceCell)      :: cell
c 
c **** declaration of internally used variables
c     nsurf,nbulk,j = loop counters
c     index         = species index
c     k             = counter
c     ier           = error flag for memory allocation
c     smodel        = copy of pointer to the surface model
      Integer :: ispec,ireact,iphase,ibulk,index,j,ier,k
      type(finiteRateSurfaceModelType),pointer :: smodel
      real*8,dimension(cell%frs_record%isrt) :: frates, brates, 
     &                 eqa, eqc, df_dtg, df_dtw, deqa_dtg, deqa_dtw,
     &                 deqc_dtg, deqc_dtw, db_dtg, db_dtw
      real*8,dimension(cell%frs_record%nspecies) :: gibbs,
     &                 dg_dtw, dg_dtg
      real*8 :: tmp, prate, rflux, pf, pb, dFdtg, dFdtw, mdotc
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- copy model pointer (easy access)
      smodel => cell%frs_record

c---- initialize data arrays
      cell%srhs      = 0.0d0
      cell%sjacobian = 0.0d0

      frates   = 0.0d0
      brates   = 0.0d0
      eqa      = 0.0d0
      eqc      = 0.0d0
      df_dtg   = 0.0d0
      df_dtw   = 0.0d0
      deqa_dtg = 0.0d0
      deqa_dtw = 0.0d0
      deqc_dtg = 0.0d0
      deqc_dtw = 0.0d0
      db_dtg   = 0.0d0
      db_dtw   = 0.0d0
      gibbs    = 0.0d0
      dg_dtw   = 0.0d0
      dg_dtg   = 0.0d0
      dFdtg    = 0.0d0
      dFdtw    = 0.0d0
      

c---- work through all forward reaction rates one at a time
      do ireact = 1,smodel%isrt

c        ---- identify gas phase participant ----      
         k = smodel%RR(ireact,1)
         if (smodel%istype(ireact) .eq. 4) k = smodel%RP(ireact,1)

c        ---- compute forward rate and derivatives ----
         call srate(cell%wall_tmpr, cell%gas_tmpr, smodel%wmass(k),
     &              smodel%sdensity(ireact), smodel%sitep(ireact),
     &              smodel%istype(ireact),
     &              smodel%Cf(ireact), smodel%eta(ireact),
     &              smodel%Eact(ireact), frates(ireact),
     &              df_dtg(ireact), df_dtw(ireact))

      enddo

c---- obtain Gibbs energies for all species either by Lewis data or backward rates
      call finite_rate_surf_gibbs_energy(smodel, gibbs, dg_dtg, dg_dtw,
     &                                   cell%gas_tmpr, cell%wall_tmpr, 
     &                                   frates, df_dtg, df_dtw)


c---- calculate activity based equilibrium constants & derivatives
      call finite_rate_surf_eqcon(smodel, gibbs, dg_dtg, dg_dtw,
     &                   cell%wall_tmpr, eqc, deqc_dtg, deqc_dtw)

c---- calculate backward rates & derivatives
      do ireact = 1,smodel%isrt
         brates(ireact) = frates(ireact)/eqc(ireact)
         db_dtg(ireact) = (df_dtg(ireact) - brates(ireact)*
     &                                    deqc_dtg(ireact))/eqc(ireact)
         db_dtw(ireact) = (df_dtw(ireact) - brates(ireact)*
     &                                    deqc_dtw(ireact))/eqc(ireact)
      enddo


c---- calculate flux for each reaction and track derivatives
      do ireact = 1,smodel%isrt

c        ----- sum forward and backward concentration product term -----
         pf = 1.0d0
         pb = 1.0d0
         do ispec = 1,smodel%nspecies
            pf = pf*cell%concentration(ispec)**smodel%vr(ireact,ispec)
            pb = pb*cell%concentration(ispec)**smodel%vp(ireact,ispec)
         enddo

c        ---- account for rates turned on and compute reaction flux ----
         pf = pf * smodel%isrfon(ireact) * smodel%sfraction(ireact)
         pb = pb * smodel%isrbon(ireact) * smodel%sfraction(ireact)
         rflux = frates(ireact)*pf - brates(ireact)*pb

c        ---- compute this contribution to species production and all derivatives ----         
         do ispec = 1,smodel%nspecies
            prate = smodel%v(ireact,ispec) * rflux

c           ---- explicit contribution ----           
            cell%srhs(ispec) = cell%srhs(ispec) + prate

c           ---- derivatives w.r.t. species concentrations ----            
c                this causes problems if concentration is exactly zero!
            do j = 1,smodel%nspecies
               cell%sjacobian(ispec,j) =  cell%sjacobian(ispec,j) + 
     &                       smodel%v(ireact,ispec) * (
     &                           frates(ireact)*pf*smodel%vr(ireact,j)
     &                         - brates(ireact)*pb*smodel%vp(ireact,j)
     &                       ) / cell%concentration(j)
            enddo

c           ---- derivatives w.r.t. temperatures ----
            dFdtg = smodel%v(ireact,ispec) * 
     &                         (pf*df_dtg(ireact)-pb*db_dtg(ireact))
            dFdtw = smodel%v(ireact,ispec) * 
     &                         (pf*df_dtw(ireact)-pb*db_dtw(ireact))

            cell%sjacobian(ispec,smodel%nspecies+1) = 
     &                  cell%sjacobian(ispec,smodel%nspecies+1) + dFdtg
            cell%sjacobian(ispec,smodel%nspecies+2) = 
     &                  cell%sjacobian(ispec,smodel%nspecies+2) + dFdtw

         enddo
         
      enddo


c---- add pyrolysis blowing term in
      if (smodel%nblwflag .eq. 1) then

         if (smodel%nebc .ge. 1) then
c           ---- sum up bulk blowing rates to add to the pyrolysis ----
c              - remember that mdot_c is the NEGATIVE of bulk      -
c              - production                                        -
            do ibulk = smodel%ibstart,smodel%ibend
               do ispec = 1,smodel%ngps
                  cell%srhs(ispec) = cell%srhs(ispec) - 
     &                              smodel%pyroadd(ispec)*
     &                              smodel%wmass(ibulk)*cell%srhs(ibulk)
               
                  do j = 1,smodel%nspecies+2  !include temperature dependencies
                     cell%sjacobian(ispec,j) = cell%sjacobian(ispec,j) -
     &                     smodel%pyroadd(ispec)*
     &                     smodel%wmass(ibulk)*cell%sjacobian(ibulk,j)
                  enddo
            
               enddo
            enddo

         else

c           ---- straight explicit mass flow addition  ----
            do ispec = 1,smodel%ngps
               cell%srhs(ispec) = cell%srhs(ispec) + 
     &                                         smodel%pyroadd(ispec)
            enddo
         
         endif

      endif

c---- diagonalize temperature rows
      cell%sjacobian(smodel%nspecies+1,smodel%nspecies+1) = 1.0d0
      cell%sjacobian(smodel%nspecies+2,smodel%nspecies+2) = 1.0d0

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine finite_rate_compute_system
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ------------Subroutine finite_rate_surf_qss_matrix()-----------------
c ---------------------------------------------------------------------
      Subroutine finite_rate_surf_qss_matrix(cell, qsslhs, 
     &                                       qssrhs, neqns)
c
c *********************************************************************
c ******    fill in the LHS and RHS terms for the QSS            ******
c ******    (steady-state) surface coverage solution.            ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      February 12, 2011                        ******
c ****** last updated:  February 12, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     cell   = cell with all concentrations already set
c     qsslhs = left-hand side Jacobian matrix
c     qssrhs = right-hand side explicit matrix
c     neqns  = number of equations (number of surface vars + number of phases)
      Integer                       :: neqns
      type(finiteRateSurfaceCell)   :: cell
      Real*8,dimension(neqns)       :: qssrhs
      Real*8,dimension(neqns,neqns) :: qsslhs
c 
c **** declaration of internally used variables
c     nsurf,nsa,j = loop counters
c     index       = species index
c     k           = counter
c     ier         = error flag for memory allocation
c     smodel      = copy of pointer to the surface model
      Integer :: nsurf,nsa,irow,icol,j,k
      type(finiteRateSurfaceModelType),pointer :: smodel
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- copy model pointer (easy access)
      smodel => cell%frs_record

c---- compute rates for the current system concentrations
      call finite_rate_compute_system(cell)

c---- clear matrices
      qsslhs = 0.0d0
      qssrhs = 0.0d0

c---- fill in LHS and RHS: surface species equations
      do k = 1,cell%num_surf
         qssrhs(k) = -cell%srhs(cell%num_gas+k)
         do j = 1,cell%num_surf
            qsslhs(k,j) = cell%sjacobian(cell%num_gas+k,cell%num_gas+j)
         enddo
      enddo

c---- fill in LHS and RHS: site conservation equations
      irow = cell%num_surf
      icol = 0
      do nsurf = 1,smodel%nsp
        do nsa = 1,smodel%nspas(nsurf)
          irow = irow + 1
          qssrhs(irow) = smodel%sdenas(nsurf,nsa)
          do k = 1,smodel%nspass(nsurf,nsa)
            icol = icol + 1
            qsslhs(irow,icol) = 1.0d0
            qsslhs(icol,irow) = 1.0d0
            qssrhs(irow) = qssrhs(irow) - 
     &                     cell%concentration(cell%num_gas+icol)
          enddo
        enddo
      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine finite_rate_surf_qss_matrix
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c -------------Function finite_rate_surf_initialize()------------------
c ---------------------------------------------------------------------
      Integer Function finite_rate_surf_initialize(cell, smodel)
c
c *********************************************************************
c ******    initialize a cell with a given model and set         ******
c ******    initial conditions for the concentrations.           ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      April    30, 2009                        ******
c ****** last updated:  June     16, 2010                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     cell   = existing cell to create and set arrays for
c     smodel = existing model to initialize from
      type(finiteRateSurfaceCell)      :: cell
      type(finiteRateSurfaceModelType),target :: smodel
c 
c **** declaration of internally used variables
c     nsurf,nbulk,nsa,j = loop counters
c     index             = species index
c     k                 = counter
c     ier               = error flag for memory allocation
      Integer :: nsurf,nbulk,nsa,index,j,ier,k
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN FUNCTION BODY ********************
c 

c---- copy array sizing (size for num_species + 2 temperatures)
      cell%num_species = smodel%nspecies
      cell%num_gas     = smodel%ngps
      cell%num_nongas  = smodel%nspecies - smodel%ngps
      cell%num_nonbulk = smodel%nspecies
      do nbulk = 1,smodel%nbp
           cell%num_nonbulk = cell%num_nonbulk - smodel%nbps(nbulk)
      enddo
      cell%num_surf = 0
      do nsurf = 1,smodel%nsp
           cell%num_surf = cell%num_surf + smodel%nsps(nsurf)
      enddo
      cell%num_eqns    = cell%num_species + 2

c---- allocate arrays for this cell
      allocate(cell%concentration(cell%num_species),
     &         cell%delta_conc(cell%num_species),
     &         cell%srhs(cell%num_eqns),
     &         cell%sjacobian(cell%num_eqns,cell%num_eqns),
     &         stat=ier)

c---- initialize all computation arrays to zero
      cell%delta_conc = 0.0d0
      cell%srhs       = 0.0d0
      cell%sjacobian  = 0.0d0

c---- store a copy to this model (pointer assignment)
      cell%frs_record => smodel

c---- concentrations for all species set to zero except surface sites
c     gas phase species will be overridden in the computation
      cell%concentration = 1.0d-15

c---- set surface densities for all surface site species (empty)      
      k=smodel%ngps
      do nsurf=1,smodel%nsp
        do nsa=1,smodel%nspas(nsurf)
          do j=1,smodel%nspass(nsurf,nsa)
            k=k+1
            if (j.eq.1) cell%concentration(k) = smodel%sdenas(nsurf,nsa)
          enddo
        enddo
      enddo

c---- set mole fractions in all bulk phases
      do nbulk=1,smodel%nbp
        do j=1,smodel%nbps(nbulk)
          k=k+1
          cell%concentration(k) = smodel%bmf(nbulk,j)
        enddo
      enddo     

      finite_rate_surf_initialize = ier
c 
c **************** END FUNCTION BODY *********************
c
      End Function finite_rate_surf_initialize
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ------------Subroutine read_finite_rate_surf_file()------------------
c ---------------------------------------------------------------------
      Subroutine read_finite_rate_surf_file(fname, smodel, ierr)
c
c *********************************************************************
c ******    read in a finite rate surface file and store the     ******
c ******    data in an object                                    ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      April     24, 2009                       ******
c ****** last updated:  May       17, 2011                       ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     fname    = name of file to open
c     smodel   = a finite rate model data record to fill in
c     ierr     = error flag for return
      Character*128                    :: fname
      type(finiteRateSurfaceModelType) :: smodel
      Integer                          :: ierr
c 
c **** declaration of internally used variables
c      i,j,k       = loop counters
c      nbulk,nsurf = loop counters for bulk and surface indices
c      nsa         = loop counter for active site index
c      nphase      = running count of number of phases
c      ier         = internal error for allocations
c      kim,kmol    = indices for activation energy
c      kf          = phase index
c      ncount      = species counter for surface species
c      Edsum       = reactant energy summation used to get Eact
      Integer :: i,j,k,nsurf,nbulk,nsa,nphase,ier,kim,kmol,kf,ncount
      Real*8  :: Edsum
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      ierr = 0
      open(10, file=fname, status='old', iostat=ierr)
	if (ierr .ne. 0) return
      rewind(10)

c-------------------------------------------------
c-----Read in all phase and species information
c-------------------------------------------------
      read(10,'(a200)') smodel%title 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) smodel%nsp,smodel%nbp
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) smodel%nblwflag,smodel%initsurf
      read(10,*) 
c     ----- declare & initialize arrays -----
      smodel%nphases = 1 + smodel%nsp + smodel%nbp
      allocate(smodel%nsps(smodel%nsp),
     &         smodel%nspas(smodel%nsp),
     &         smodel%sfrc(smodel%nsp),
     &         smodel%sden(smodel%nsp),
     &         smodel%ithermo(smodel%nsp),
     &         smodel%den(smodel%nbp),
     &         smodel%nbps(smodel%nbp),
     &         smodel%por(smodel%nbp),
     &         smodel%vfb(smodel%nbp),
     &         smodel%phases(smodel%nphases),stat=ier)
c     ----- read Gas phase declaration -----
      nphase=1
      read(10,*) 
      read(10,*)  
      read(10,*) smodel%phases(nphase), smodel%ngps
      read(10,*)
      smodel%nspecies = smodel%ngps 
c     ----- read Surface phase declarations -----
      read(10,*) 
      read(10,*) 
      read(10,*)
      smodel%nspasmax=0
      do nsurf=1,smodel%nsp
         nphase=nphase+1
         read(10,*) smodel%phases(nphase),smodel%sfrc(nsurf),
     &              smodel%nspas(nsurf),smodel%ithermo(nsurf)
         if ( smodel%nspas(nsurf) .gt. smodel%nspasmax) 
     &                    smodel%nspasmax = smodel%nspas(nsurf)
      enddo
      read(10,*) 
c     ----- declare & initialize nspass & sdenas arrays -----
      allocate(smodel%nspass(smodel%nsp,smodel%nspasmax),
     &         smodel%sdenas(smodel%nsp,smodel%nspasmax),stat=ier)
      smodel%nspass  = 0
      smodel%sdenas  = 0.0d0
      smodel%sden    = 0.0d0
      smodel%nsps    = 0
      smodel%sdentot = 0.0d0
      smodel%sfrcI   = 0.0d0
c     ----- read Surface Active Site declarations -----
      read(10,*) 
      read(10,*) 
      read(10,*)
      do nsurf=1,smodel%nsp
         do nsa = 1,smodel%nspas(nsurf)
           read(10,*) smodel%sdenas(nsurf,nsa),
     &                smodel%nspass(nsurf,nsa)
           smodel%nsps(nsurf) = smodel%nsps(nsurf) + 
     &                          smodel%nspass(nsurf,nsa)
           smodel%sden(nsurf) = smodel%sden(nsurf) + 
     &                          smodel%sdenas(nsurf,nsa)
         enddo
         smodel%nspecies = smodel%nspecies + smodel%nsps(nsurf)
         smodel%sdentot = smodel%sdentot + 
     &                    smodel%sfrc(nsurf)*smodel%sden(nsurf)
         if (smodel%nspas(nsurf) .eq. 0)
     &         smodel%sfrcI = smodel%sfrcI + smodel%sfrc(nsurf)
      enddo
      read(10,*)
c     ----- read Bulk phase declarations -----
      smodel%vfbI    = 0.0d0
      smodel%bdentot = 0.0d0
      smodel%bportot = 0.0d0
      if (smodel%nbp .gt. 0) then
         read(10,*) 
         read(10,*) 
         read(10,*) 
         do nbulk=1,smodel%nbp
           nphase=nphase+1
           read(10,*) smodel%phases(nphase),smodel%den(nbulk),
     &                smodel%por(nbulk),smodel%vfb(nbulk),
     &                smodel%nbps(nbulk)
           smodel%nspecies = smodel%nspecies + smodel%nbps(nbulk)
           if (smodel%nbps(nbulk) .eq. 0)
     &           smodel%vfbI = smodel%vfbI + smodel%vfb(nbulk)
           smodel%bdentot = smodel%bdentot + 
     &                      smodel%vfb(nbulk)*smodel%den(nbulk)
           smodel%bportot = smodel%bportot + 
     &                      smodel%vfb(nbulk)*smodel%por(nbulk)
         enddo
         read(10,*) 
      endif
c     ----- compute start/end indices -----
      smodel%igstart = 1
      smodel%igend   = smodel%ngps 
      smodel%isstart = smodel%igend + 1
      ncount = 0
      do nsurf=1,smodel%nsp
         ncount = ncount + smodel%nsps(nsurf)
      enddo
      smodel%isend   = smodel%isstart + ncount - 1
      smodel%ibstart = smodel%isstart + ncount
      smodel%ibend   = smodel%nspecies

c-----------------------------------------------
c-----Read in species lists       --------------
c-----------------------------------------------
c     ----- declare & initialize arrays -----
      allocate(smodel%species(smodel%nspecies),
     &         smodel%wmass(smodel%nspecies),
     &         smodel%Ediss(smodel%nspecies),
     &         smodel%Ed(smodel%nspecies),
     &         smodel%kphase(smodel%nspecies),
     &         smodel%kthermo(smodel%nspecies),
     &         smodel%ksites(smodel%nspecies),
     &         smodel%bmf(smodel%nbp,smodel%nspecies),
     &         stat=ier)
      smodel%kphase  = 0
      smodel%ksites  = 0
      smodel%kthermo = 1
      smodel%Ediss   = 0.0d0
      smodel%Ed      = 0.0d0
      smodel%bmf     = 0.0d0
c     ----- read Gas phase species list -----
      read(10,*) 
      read(10,*) 
      do k=1,smodel%ngps
        read(10,*) smodel%species(k),smodel%wmass(k),smodel%Ediss(k)
        smodel%kphase(k)=1
      enddo
      k = smodel%ngps
      read(10,*) 
c     ----- read Surface phase species list -----
      read(10,*) 
      read(10,*)
      do nsurf=1,smodel%nsp
        do nsa = 1,smodel%nspas(nsurf)
          do j=1,smodel%nspass(nsurf,nsa)
            k=k+1
            read(10,*) smodel%species(k),smodel%wmass(k),smodel%Ed(k)
            smodel%kphase(k)=1+nsurf
            smodel%ksites(k)=nsa
            smodel%kthermo(k)=smodel%ithermo(nsurf)
          enddo
        enddo
      enddo
      read(10,*)
c     ----- read Bulk phase species list -----
      if (smodel%nbp .gt. 0) then
        read(10,*)
        read(10,*)
        do nbulk=1,smodel%nbp 
          do j=1,smodel%nbps(nbulk)
            k=k+1
            read(10,*) smodel%species(k),smodel%wmass(k),
     &                 smodel%bmf(nbulk,j)
            smodel%kphase(k)=1+smodel%nsp+nbulk
          enddo
        enddo
        read(10,*)
      endif

c-----------------------------------------------
c-----Read in reaction information--------------
c-----------------------------------------------
      read(10,*) 
      read(10,*)
      read(10,*) smodel%isrt
      read(10,*) 
c     ----- declare & initialize arrays -----
c           ** vr,vp,v are indexed from 0 as a catch for no species
      allocate(smodel%RR(smodel%isrt,ISMODEL_MAX_REACT),
     &         smodel%RP(smodel%isrt,ISMODEL_MAX_PROD),
     &         smodel%RVR(smodel%isrt,ISMODEL_MAX_REACT),
     &         smodel%RVP(smodel%isrt,ISMODEL_MAX_PROD),
     &         smodel%vr(smodel%isrt,0:smodel%nspecies),
     &         smodel%vp(smodel%isrt,0:smodel%nspecies),
     &         smodel%v(smodel%isrt,0:smodel%nspecies),stat=ier)
      smodel%RR  = 0
      smodel%RP  = 0
      smodel%RVR = 0
      smodel%RVP = 0
      smodel%vp  = 0
      smodel%vr  = 0
      smodel%v   = 0
c---- Reactant and product species
      read(10,*) 
      do i=1,smodel%isrt
        read(10,*) smodel%RR(i,1),smodel%RR(i,2),smodel%RR(i,3),
     &             smodel%RP(i,1),smodel%RP(i,2),smodel%RP(i,3)
      enddo
      read(10,*) 
c---- Reactant and product stoichiometric coefficients
      read(10,*) 
      do i=1,smodel%isrt
        read(10,*) smodel%RVR(i,1),smodel%RVR(i,2),smodel%RVR(i,3),
     &             smodel%RVP(i,1),smodel%RVP(i,2),smodel%RVP(i,3)
        do j=1,ISMODEL_MAX_REACT
          smodel%vr(i,smodel%RR(i,j)) = smodel%vr(i,smodel%RR(i,j))
     &                                + smodel%RVR(i,j)
        enddo
        do j=1,ISMODEL_MAX_PROD  
          smodel%vp(i,smodel%RP(i,j)) = smodel%vp(i,smodel%RP(i,j))
     &                                + smodel%RVP(i,j)
        enddo
        do k=1,smodel%nspecies
          smodel%v(i,k)=smodel%vp(i,k)-smodel%vr(i,k)
        enddo
      enddo
      read(10,*) 

c-----------------------------------------------------
c---- Read in rate constant information --------------
c-----------------------------------------------------
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*)
	read(10,*)
	read(10,*) 
c     ----- declare & initialize arrays -----
      allocate(smodel%istype(smodel%isrt),
     &         smodel%Cf(smodel%isrt),
     &         smodel%eta(smodel%isrt),
     &         smodel%Ea(smodel%isrt),
     &         smodel%isrfon(smodel%isrt),
     &         smodel%isrbon(smodel%isrt),stat=ier)
      smodel%istype = 0
      smodel%Cf     = 0.0d0
      smodel%eta    = 0.0d0
      smodel%Ea     = 0.0d0
      smodel%isrfon = 1
      smodel%isrbon = 1
c     ----- read parameters -----
      do i=1,smodel%isrt
        read(10,*) smodel%istype(i),smodel%Cf(i),smodel%eta(i),
     &             smodel%Ea(i),smodel%isrfon(i),smodel%isrbon(i)
      enddo
      read(10,*) 

c--------------------------------------------------------------
c---- For surface phases without thermochemical data ----------
c---- Read in desorption reactions or equilibrium constants ---
c--------------------------------------------------------------
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
      read(10,*) 
c     ----- declare & initialize arrays -----
      allocate(smodel%iadtype(smodel%isrt),
     &         smodel%iadform(smodel%isrt),
     &         smodel%Cfad(smodel%isrt),
     &         smodel%etaad(smodel%isrt),
     &         smodel%vad(smodel%isrt),
     &         smodel%Edes(smodel%isrt),
     &         smodel%kadsr(smodel%nspecies),stat=ier)
      smodel%iadtype = 0
      smodel%iadform = 0
      smodel%kadsr   = 0
      smodel%Cfad    = 0.0d0
      smodel%etaad   = 0.0d0
      smodel%etaad   = 0.0d0
      smodel%Edes    = 0.0d0
      smodel%vad     = 0.0d0
      do i=1,smodel%isrt
c------ If an adsorption reaction
        if ((smodel%istype(i).eq.1).or.(smodel%istype(i).eq.5)) then
c-------- identify surface species
          k=smodel%RP(i,1)
          smodel%kadsr(k)=i
c-------- If no thermodynamic data for surface species k
          if(smodel%kthermo(k).eq.0) then
            read(10,*) smodel%iadtype(i),smodel%iadform(i),
     &                 smodel%Cfad(i),   smodel%etaad(i),
     &                 smodel%vad(i),    smodel%Edes(i)
            smodel%Ed(k)=smodel%Edes(i)
          endif
        endif
      enddo

c---- compute an activation energy & total site density for all reactions
      allocate(smodel%Eact(smodel%isrt),
     &         smodel%sitep(smodel%isrt),
     &         smodel%sdensity(smodel%isrt),
     &         smodel%sfraction(smodel%isrt),
     &         smodel%pyroadd(smodel%ngps),stat=ier)
      smodel%Eact      = 0.0d0
      smodel%sitep     = 0.0d0
      smodel%sdensity  = 0.0d0
      smodel%sfraction = 0.0d0
      smodel%pyroadd   = 0.0d0

      do i=1,smodel%isrt
        smodel%Eact(i)=smodel%Ea(i)
        if (smodel%istype(i).eq.3) then

          Edsum = 0.0d0
          do j=1,ISMODEL_MAX_REACT
             k=smodel%RR(i,j)
             if (k .gt. 0)
     &          Edsum = Edsum + smodel%RVR(i,j)*smodel%Ed(k)
          enddo

          kmol=smodel%RP(i,1)
          smodel%Eact(i) = max(smodel%Ea(i),
     &                Edsum-smodel%Ediss(kmol))
        endif

        do j=1,ISMODEL_MAX_REACT
           k=smodel%RR(i,j)
	     if (k .gt. 0) then
             kf = smodel%kphase(k)
             if ( (kf .gt. 1) .and. (kf .le. (smodel%nsp+1)) ) then
               smodel%sitep(i)     = smodel%sitep(i) 
     &                               + dble(smodel%RVR(i,j))
               smodel%sdensity(i)  = smodel%sden(kf-1)
               smodel%sfraction(i) = smodel%sfrc(kf-1)
             endif
           endif
        enddo

      enddo

c---- compute stochiometric coefficient sum for all reactions
      allocate(smodel%sumv(smodel%isrt,smodel%nphases),stat=ier)
      smodel%sumv = 0.0d0
      do i=1,smodel%isrt
        do k=1,smodel%isend
           j = smodel%kphase(k)
           smodel%sumv(i,j) = smodel%sumv(i,j) + smodel%v(i,k)
        enddo
      enddo

c---- debugging off
      smodel%ifrm_debug = 0

c---- blowing off until a pyrolysis file can be read
      smodel%isblwon = 0
      smodel%nblw    = 0
      smodel%cyield  = 1.0d0
      smodel%ha      = 0.0d0
      
c---- close file and return -----
      close (10)


c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine read_finite_rate_surf_file
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ------------Subroutine read_finite_rate_blow_file()------------------
c ---------------------------------------------------------------------
      Subroutine read_finite_rate_blow_file(blwname, smodel, ierr)
c
c *********************************************************************
c ******    read in a finite rate surface blowing data file      ******
c ******    and store into an existing object                    ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      June     11, 2010                        ******
c ****** last updated:  April    19, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     blwname  = name of file to open
c     smodel   = a finite rate model data record to fill in
c     ierr     = error flag for return
      Character*128                    :: blwname
      type(finiteRateSurfaceModelType) :: smodel
      Integer                          :: ierr
c 
c **** declaration of internally used variables
c      n,k  = loop counters
c      cval = value of [1/cyield - 1]
c      sumxkk = sum of mole fraction x molecular weight
      Integer :: n,k
      real*8  :: cval, sumxkk
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      ierr = 0

c---- check if we need to do this ----
      if (smodel%nblwflag .le. 0) return

c---- open file ----
      open(10, file=blwname, status='old', iostat=ierr)
      if (ierr .ne. 0) return
      rewind(10)

c---- read the header of the blowing file ----
      read(10,*)
      read(10,*)
      read(10,*) smodel%nebc
      if ((smodel%nebc .lt. 0) .or. (smodel%nebc .gt. 2)) then
         ierr = 1
         return
      endif
      if (smodel%nebc .eq. 1) then
         read(10,*) smodel%cyield
         smodel%nblw = 1
         smodel%ha   = 0.0d0
      elseif (smodel%nebc .eq. 2) then
         read(10,*) smodel%cyield, smodel%ha
         smodel%nblw = 1
      else
         read(10,*) smodel%nblw
         smodel%cyield = 1.0d0    !cyield must be 1.0 to not affect bulk production
         smodel%ha     = 0.0d0
      endif

c---- allocate memory ----
      allocate(smodel%rmdotg(smodel%nblw),
     &         smodel%blwmf(smodel%nblw,smodel%ngps), stat=ierr)
      smodel%rmdotg  = 0.0d0
      smodel%blwmf   = 0.0d0

c---- read the body of the blowing file ----
      if (smodel%nebc .eq. 0)
     &   read(10,*) (smodel%rmdotg(n), n=1,smodel%nblw)

      do k=1,smodel%ngps
         read(10,*) (smodel%blwmf(n,k), n=1,smodel%nblw)
      enddo

c---- compute the effective mass addition term (save operations) ----
      if (smodel%nebc .gt. 0) then
         cval = 0.0d0
         if (smodel%cyield .gt. 0.0d0)
     &           cval = 1.0d0/smodel%cyield - 1.0d0
         sumxkk = 0.0d0
         do k=1,smodel%ngps
            sumxkk = sumxkk + smodel%blwmf(1,k)*smodel%wmass(k)
         enddo
         if (sumxkk .eq. 0.0d0) sumxkk = 1.0d0
         do k=1,smodel%ngps
            smodel%pyroadd(k) = smodel%blwmf(1,k)*cval/sumxkk
         enddo
         smodel%rmdotg(1) = 0.0d0  !explicit mass addition
      else
         do n=1,smodel%nblw
            sumxkk = 0.0d0
            do k=1,smodel%ngps
               sumxkk = sumxkk + smodel%blwmf(n,k)*smodel%wmass(k)
            enddo
            if (sumxkk .eq. 0.0d0) sumxkk = 1.0d0
            do k=1,smodel%ngps
               smodel%pyroadd(k) = smodel%pyroadd(k) + 
     &                       smodel%blwmf(n,k)*smodel%rmdotg(n)/sumxkk
            enddo
         enddo
      endif

c---- flag blowing on ----
      smodel%isblwon = 1

c---- close file and return ----
      close (10)

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine read_finite_rate_blow_file
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c --------------Subroutine write_finite_rate_surf()--------------------
c ---------------------------------------------------------------------
      Subroutine write_finite_rate_surf(smodel, iunit, ierr)
c
c *********************************************************************
c ******    write a model back to STDOUT as a diagnostic tool    ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      November 12, 2009                        ******
c ****** last updated:  January  31, 2013                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     smodel   = a finite rate model to write
c     iunit    = record unit to write to (must be open)
c     ierr     = error flag for return
      type(finiteRateSurfaceModelType) :: smodel
      Integer                          :: iunit,ierr
c 
c **** declaration of internally used variables
c      j,k,m,n  = loop counters
c      Edes     = desorption energy value
c      sline    = output string line
      Integer       :: j,k,m,n
	Real*8        :: Edes
	character*200 :: sline
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      ierr = 0

      write(iunit,*) trim(smodel%title)
      write(iunit,*) '------------------------------------------------'
      write(iunit,*) 'Number of Surface and bulk phases',
     &               ' (gas phase=1 by definition)'
      write(iunit,*) 'nsp, nbp'
      write(iunit,'(i4,x,i4)') smodel%nsp, smodel%nbp
      write(iunit,*)

      write(iunit,'(a,a)') ' Blowing/pyrolyzing gas flows? (0=NO, ',
     &       '1=Yes) and Surface Initialization (0=empty sites, 1=QSS)'
      write(iunit,*) 'nblwflag     initsurf'
      write(iunit,'(x,i1,12x,i1)') smodel%nblwflag,smodel%initsurf
      write(iunit,*)

      write(iunit,*) 'Number of gas phase species participating',
     &               ' in surface reactions'
      write(iunit,*) 'Name      ngps  #Phase'
      write(iunit,'(x,a,x,i3,4x,i3)')
     &                       smodel%phases(1), smodel%ngps, 1
      write(iunit,*)
	
      write(iunit,*) 'For each surface phase: list name, ',
     &               'surface fraction and'
      write(iunit,*) 'number of active site sets, and ',
     &               'thermo availability (0=No,1=Yes)'
      write(iunit,*) 'Name       sfrc     nspas  iThermo  #Phase'
      do n = 1,smodel%nsp
         write(iunit,'(x,a,x,f5.3,2x,i3,6x,i1,6x,i3)') 
     &                smodel%phases(n+1),smodel%sfrc(n),
     &                smodel%nspas(n),smodel%ithermo(n),n+1
	enddo
      write(iunit,*)

      write(iunit,*) 'For each surface phase with 1 or more sets',
     &               ' of active sites, list'
      write(iunit,*) 'the site density and number of',
     &               'species for each active site set'
      write(iunit,*) 'sdenas (mol/m2)      nspass     phase#/site#'
      do n=1,smodel%nsp
        do k=1,smodel%nspas(n)
          write(iunit,'(x,es12.5,7x,i3,10x,i1,a,i3)') 
     &                smodel%sdenas(n,k),
     &                smodel%nspass(n,k),n+1,'/',k
        enddo
      enddo
      write(iunit,*)

      if (smodel%nbp .gt. 0) then
         write(iunit,'(a,a)') ' For each bulk phase: list name, ',
     &                  'density, porosity, '
         write(iunit,*) 'volume fraction and number of bulk species '
         write(iunit,'(x,a,a)') 'Name       sden         porosity',
     &                  '  vol. fract.  nbps   #Phase'
         do n = 1,smodel%nbp
            write(iunit,'(x,a,x,es10.4,2x,f7.4,3x,f7.4,5x,i3,4x,i3)') 
     &                     smodel%phases(smodel%nsp+1+n),
     &                     smodel%den(n),smodel%por(n),smodel%vfb(n),
     &                     smodel%nbps(n),n+1+smodel%nsp
	   enddo
         write(iunit,*)
      endif
	
      write(iunit,*) 'Order of gas species'
      write(iunit,*) 'Name        Molar mass     Ediss          ',
     &               '#Species'
      do n = smodel%igstart,smodel%igend
         write(iunit,'(x,a,x,f11.8,4x,es11.4,4x,i3)')
     &      smodel%species(n),smodel%wmass(n),smodel%Ediss(n),n
      enddo
      write(iunit,*)
      k = smodel%igend + 1

      write(iunit,'(x,a,a)') 'Order of surface species for ',
     &                       ' (number each consecutively)'
      write(iunit,'(x,a,a)') 'Name   Molar mass     Edes           ',
     &                       'Species#   Phase#/Site#'
      do j = 1,smodel%nsp
         do m = 1,smodel%nspas(j)
            do n = 1,smodel%nspass(j,m)
               Edes = 0.0d0
               if (smodel%kadsr(k) .gt. 0)
     &                       Edes = smodel%Edes( smodel%kadsr(k) )
               write(iunit,'(x,a,x,f11.8,4x,es11.4,4x,i3,7x,i3,a,i3)')
     &                       trim(smodel%species(k)),smodel%wmass(k),
     &                       Edes,k,j,'/',m
               k = k + 1
            enddo
         enddo
      enddo
      write(iunit,*)

      if (smodel%nbp .gt. 0) then
         write(iunit,'(x,a,a)') 'Order of bulk species ',
     &                  '(number each consecutively)'
         write(iunit,'(x,a,a)') 'Name   Molar mass     bmf        ',
     &                       '    Species#   Phase#'
         do j = 1,smodel%nbp
           do n = 1,smodel%nbps(j)
             write(iunit,'(x,a,x,f11.8,4x,es11.4,4x,i3,7x,i3)')
     &                   trim(smodel%species(k)),smodel%wmass(k),
     &                   smodel%bmf(j,n),k,j
             k = k + 1
           enddo
         enddo
         write(iunit,*)
	endif

      write(iunit,*) 'Total number of surface reactions'
      write(iunit,*) 'nsrt'
      write(iunit,'(i4)') smodel%isrt
      write(iunit,*) 

      write(iunit,*) 'Reactant/product species for each ',
     &               'forward surface reaction'
      do n = 1,smodel%isrt
         write(iunit,101) smodel%RR(n,:),smodel%RP(n,:),n
      enddo
      write(iunit,*) 

      write(iunit,*) 'Stoichiometric coefficients for each ',
     &               'surface reaction'
      do n = 1,smodel%isrt
         write(iunit,101) smodel%RVR(n,:),smodel%RVP(n,:),n
      enddo
      write(iunit,*) 

      write(iunit,*) 'Reaction parameters for each type of reaction:'
      write(iunit,*) 
     &     'Arrhenius:             0,   Cf, beta,   Ea, isrfon, isrbon'
      write(iunit,*) 
     &     'Adsorption:            1,   S0, beta, Eads, isrfon, isrbon'
      write(iunit,*) 
     &     'Eley-Rideal:           2,  Ger, beta,  Eer, isrfon, isrbon'
      write(iunit,*) 
     &     'Langmuir-Hinschelwood: 3,  Clh, beta,   Em, isrfon, isrbon'
      write(iunit,*) 
     &     'Sublimation:           4, a0Pv, beta, Esub, isrfon, isrbon'
      write(iunit,*) 
     &     'Arrhenius Adsorption:  5,   Cf, beta,   Ea, isrfon, isrbon'
	write(iunit,*) 
     &     'Type, (Param), beta,   Ea, isrfon, isrbon'
	do n=1,smodel%isrt
         write(iunit,102) 
     &                smodel%istype(n),smodel%Cf(n),smodel%eta(n),
     &                smodel%Ea(n),smodel%isrfon(n),smodel%isrbon(n), n
      enddo
      write(iunit,*)

      write(iunit,*) 'Desorption reaction or equilibrium',
     &               ' constant parameters:'
      write(iunit,*) 'Type 1: Desorption:'
      write(iunit,*) '        Form 0: Arrhenius'
      write(iunit,*) '        Form 1: Constant attempt frequency'
      write(iunit,*) '        Form 2: Simple transition state theory'
      write(iunit,*) '        Form 3: Complex transition state theory'
      write(iunit,*) 'Type 2: Equilibrium:'
      write(iunit,*) '        Form 0: Arrhenius'
      write(iunit,*) '        Form 1: Immobile adsorption -',
     &               ' simple transition state theory'
      write(iunit,*) '        Form 2: Immobile adsorption -',
     &               ' complex transition state theory'
      write(iunit,*) '        Form 3: Mobile adsorption -',
     &               ' simple transition state theory'
      write(iunit,*) '        Form 4: Mobile adsorption -',
     &               ' complex transition state theory'
      write(iunit,*) 'Type, Form, Cf, eta, vdes, Edes'
      do n=1,smodel%isrt
         if ((smodel%istype(n).eq.1).or.(smodel%istype(n).eq.5)) then
             write(iunit,103) 
     &                 smodel%iadtype(n),smodel%iadform(n),
     &                 smodel%Cfad(n),   smodel%etaad(n),
     &                 smodel%vad(n),    smodel%Edes(n), n
         endif
      enddo
      write(iunit,*) 

c     ----- delimiter -----
      write(iunit,*) '     --------------------     '

c     ----- not in input file, but also write inert surface data -----
      write(iunit,*)
      write(iunit,'(a,i1,a,i1)') 'Finite Rate Model Version ',
     &            MAJOR_FRM_VERSION,'-',MINOR_FRM_VERSION
      write(iunit,'(a,f7.4)') 'Inert Surface Fraction   ',smodel%sfrcI
      write(iunit,'(a,es10.4)') 'Surface Density           ',
     &                          smodel%sdentot
      write(iunit,'(a,f7.4)') 'Inert Bulk Fraction      ',smodel%vfbI
      write(iunit,'(a,es10.4)') 'Bulk density              ',
     &                          smodel%bdentot
      write(iunit,'(a,f7.4)') 'Bulk Porosity            ',smodel%bportot
      write(iunit,*)

c     ----- not in input file, but also write blowing data -----
      if (smodel%isblwon .gt. 0) then
        write(iunit,*) '     --------------------     '
        write(iunit,*)
        write(iunit,'(a,a)') ' Blowing gas Data for ',trim(smodel%title)
        write(iunit,*) '----------------------------------------------'
        write(iunit,'(x,i1,6x,a,a)') smodel%nebc,
     &            '| EBC type (0: mdotg (kg/s)) (1: cyield)',
     &            ' (2: cyield, ha (J/kg))'
        if (smodel%nebc .gt. 0) then
          write(iunit,'(x,f7.4,2x,es12.5,6x,a)') smodel%cyield, 
     &                smodel%ha,'| yield/enthalpy'      
        else
          write(iunit,'(x,i1,6x,a)') smodel%nblw,
     &                '| Number of blowing/pyrolyzing flows'
          write(sline,'(15es11.4)')
     &                             (smodel%rmdotg(n), n=1,smodel%nblw)
          sline = trim(sline)//'   | each mdotg (kg/s)'
	    write(iunit,'(a)') trim(sline)
        endif
        do k=1,smodel%ngps
          write(sline,'(15f7.4)') (smodel%blwmf(n,k), n=1,smodel%nblw)
          if (k .eq. 1) then
             sline = trim(sline)//'    '//trim(smodel%species(n))//
     &                   '  | <= begin Mole fractions of blowing gas'
          else
             sline = trim(sline)//'    '//trim(smodel%species(k))
          endif
          write(iunit,'(a)') trim(sline)
	  enddo
        write(iunit,*)
      endif

c     ----- delimiter -----
      write(iunit,*) '     --------------------     '

c     ----- not in input file, but also write Gibbs data -----
      write(iunit,*)
      write(iunit,*) 'Gibbs Free-energy data from Lewis database'
      do n=1,smodel%nspecies
        write(iunit,'(a,i3,i3,3es15.6)') trim(smodel%species(n)),
     &                           smodel%rlewis_nrange(n),0,
     &                           smodel%wmass(n)*1.0d3, 
     &                           0.0d0, 0.0d0
        if (smodel%kthermo(n) .gt. 0) then
           do j = 1,smodel%rlewis_nrange(n)
              write(iunit,104) smodel%rlewis_trange(1,j,n),
     &                         smodel%rlewis_trange(2,j,n)
              write(iunit,105) smodel%rlewis_coeffs(1:5,j,n)
		  write(iunit,105) smodel%rlewis_coeffs(6:10,j,n)
           enddo
        else
           write(iunit,*) '    Not using Gibbs data!'
        endif
      enddo

 101  format(i2,',',i2,',',i2,',',i2,',',i2,',',i2,6x,'#',i3)
 102  format(x,i1,', ',es11.4,', ',es11.4,', ',es11.4,', ',i1,', ',i1,
     &       '  #',i3)
 103  format(x,i1,', ',i1,', ',es11.4,', ',es11.4,', ',
     &       es11.4,', ',es11.4,'  #',i3)
 104  format('    ',es14.3,es14.3)
 105  format('    ',5es14.6)
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine write_finite_rate_surf
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ---------------Subroutine read_lewis_bulk_file()---------------------
c ---------------------------------------------------------------------
      Subroutine read_lewis_bulk_file(fname, smodel, gcof, gtmpr, 
     &                                ngfit, ierr)
c
c *********************************************************************
c ******    read in data from the bulk phase Lewis file          ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      April    24, 2009                        ******
c ****** last updated:  May       2, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     fname    = name of file to open
c     smodel   = a finite rate model data record to fill in
c     ierr     = error flag for return
c     gcof     = gas phase Gibbs coefficients (must be given)
c     gtmpr    = array of temperature breakpoints for gcof
c     ngfit    = number of fits for gas phase (usually 3 or 5)
      Character*200                          :: fname
      type(finiteRateSurfaceModelType)       :: smodel
      Integer                                :: ierr, ngfit
      real*8,dimension(10,ngfit,smodel%ngps) :: gcof
      real*8,dimension(ngfit+1)              :: gtmpr
c 
c **** declaration of internally used variables
c      i,j,m,n,k            = loop and index counters
c      ier                  = internal error for allocations
c      h,s,cp               = Lewis functions
c      t1,t2,t3,t4,ti,ts,tl = temperature functions for Lewis evaluations
c      rmwgt                = molecular weight (Not Currently Used)
c      hform298             = heat of formation at 298.15K (Not Currently Used)
c      dhform               = realative heat of formation  (Not Currently Used)
c      tmin,tmax            = temperature ranges
c      coeff                = Lewis coefficients from file
c      sname                = species name from file
c      nrange               = number of ranges from file
c      itype                = type of file data (0=bulk, >0=adsorbed surface)
c      imatch               = name match flag
c      ilp                  = position index
c      
      Integer     :: nrange,itype,imatch,ilp,i,j,m,n,k,ier
      real*8      :: rmwgt,hform298,dhform,tmin,tmax,coeff(10),
     &               t1,t2,t3,t4,ti,ts,tl,h,s,cp
      Character*8 :: sname
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      ierr = 0
      open(10, file=fname, status='old')
      rewind(10)

c-------------------------------------------------
c-----declare Lewis arrays
c-------------------------------------------------
      allocate(smodel%rlewis_coeffs(NLEWIS_COEFFS,NLEWIS_FITS,
     &                                       smodel%nspecies),
     &         smodel%rlewis_trange(2,NLEWIS_FITS,smodel%nspecies),
     &         smodel%rlewis_nrange(smodel%nspecies),
     &         stat=ier)

      smodel%rlewis_coeffs = 0.0d0
      smodel%rlewis_trange = 0.0d0
      smodel%rlewis_nrange = 0

c--------------------------------------------------------
c-----obtain gas phase Gibbs data from existing data
c--------------------------------------------------------

      do i=1,smodel%ngps
         smodel%rlewis_nrange(i) = ngfit
         do j = 1,ngfit
            smodel%rlewis_trange(1,j,i) = gtmpr(j)
            smodel%rlewis_trange(2,j,i) = gtmpr(j+1)
            do k = 1,10
               smodel%rlewis_coeffs(k,j,i) = gcof(k,j,i)
            enddo
         enddo
      enddo

c--------------------------------------------------------
c-----most surface phase species do not have Gibbs data
c-----however, set empty surface sites to use Gibbs data
c-----as this will automatically zero it out
c--------------------------------------------------------

      k = smodel%ngps
      do n=1,smodel%nsp
        do i=1,smodel%nspas(n)
          do j=1,smodel%nspass(n,i)
            k=k+1
            if (j .eq. 1) then
               smodel%kthermo(k)=1
               smodel%rlewis_trange(1,1,k) = 0.0d0
               smodel%rlewis_trange(2,1,k) = 1.0d+10
               smodel%rlewis_nrange(k)     = 1
            endif
          enddo
        enddo
      enddo

c--------------------------------------------------------
c-----pass through file looking for matching bulk species
c--------------------------------------------------------

      do k = 1,1000
      
	   read(10,*) sname,nrange,itype,rmwgt,hform298,dhform
      
         if (trim(sname) .eq. 'END') goto 121

c        ----- look for a bulk species that matches this name -----
c              strip off the (b) to check in the database
         do i = 1,nrange
            read(10,*)  tmin,tmax
            read(10,*)  coeff(1:5)
            read(10,*)  coeff(6:10)
            do n = smodel%ibstart,smodel%ibend
	         ilp = 10
	         do m=1,10
	            if (smodel%species(n)(m:m) .eq. '(') ilp = m 
	         enddo
               if ( (trim(sname) .eq. smodel%species(n)(1:ilp-1)) .and.
     &              (itype .eq. 0) .and. (i .lt. NLEWIS_FITS) ) then
                  smodel%rlewis_trange(1,i+1,n) = tmin
                  smodel%rlewis_trange(2,i+1,n) = tmax
                  do m = 1,10
                     smodel%rlewis_coeffs(m,i+1,n) = coeff(m)
                  enddo
                  smodel%rlewis_nrange(n)     = nrange + 2
               endif
            enddo
         enddo
	
	enddo

c----- close file -----
 121  close (10)


c-------------------------------------------------
c-----extrapolate non-gas species with constant Cp
c-----this has already been done for gas species
c-------------------------------------------------

      do n=smodel%ngps+1,smodel%nspecies
         if (smodel%rlewis_nrange(n) .gt. 1) then
c           ----- first range -----
            m = 2
            smodel%rlewis_trange(1,m-1,n) = 0.0d0
            smodel%rlewis_trange(2,m-1,n) = smodel%rlewis_trange(1,m,n)

            t1 = smodel%rlewis_trange(1,m,n)
            t2 = t1*t1
            t3 = t2*t1
            t4 = t3*t1
            ti = 1.0d0/t1
            ts = ti*ti
            tl = log(t1)
            h  = smodel%rlewis_coeffs(3,m,n)
     &         + smodel%rlewis_coeffs(2,m,n)*ti*tl
     &         - smodel%rlewis_coeffs(1,m,n)*ts
     &         + 0.5*smodel%rlewis_coeffs(4,m,n)*t1
     &         + 1.0d0/3.0d0*smodel%rlewis_coeffs(5,m,n)*t2
     &         + 0.25*smodel%rlewis_coeffs(6,m,n)*t3
     &         + 0.2*smodel%rlewis_coeffs(7,m,n)*t4
     &         + smodel%rlewis_coeffs(9,m,n)*ti
            s  = smodel%rlewis_coeffs(3,m,n)*tl
     &         - smodel%rlewis_coeffs(2,m,n)*ti
     &         - 0.5*smodel%rlewis_coeffs(1,m,n)*ts
     &         + smodel%rlewis_coeffs(4,m,n)*t1
     &         + 0.5*smodel%rlewis_coeffs(5,m,n)*t2
     &         + 1.0d0/3.0d0*smodel%rlewis_coeffs(6,m,n)*t3
     &         + 0.25*smodel%rlewis_coeffs(7,m,n)*t4
     &         + smodel%rlewis_coeffs(10,m,n)
            cp = smodel%rlewis_coeffs(1,m,n)*ts
     &         + smodel%rlewis_coeffs(2,m,n)*ti
     &         + smodel%rlewis_coeffs(3,m,n)
     &         + smodel%rlewis_coeffs(4,m,n)*t1
     &         + smodel%rlewis_coeffs(5,m,n)*t2
     &         + smodel%rlewis_coeffs(6,m,n)*t3
     &         + smodel%rlewis_coeffs(7,m,n)*t4

            smodel%rlewis_coeffs(3,m-1,n)  = cp
            smodel%rlewis_coeffs(9,m-1,n)  = (h-cp)*t1
            smodel%rlewis_coeffs(10,m-1,n) = s - cp*tl

c           ----- end range -----
            m = smodel%rlewis_nrange(n) - 1
            if (m .ge. NLEWIS_FITS) m = NLEWIS_FITS - 1
            smodel%rlewis_trange(1,m+1,n) = smodel%rlewis_trange(2,m,n)
            smodel%rlewis_trange(2,m+1,n) = 1.0d+20

            t1 = smodel%rlewis_trange(2,m,n)
            t2 = t1*t1
            t3 = t2*t1
            t4 = t3*t1
            ti = 1.0d0/t1
            ts = ti*ti
            tl = log(t1)
            h  = smodel%rlewis_coeffs(3,m,n)
     &         + smodel%rlewis_coeffs(2,m,n)*ti*tl
     &         - smodel%rlewis_coeffs(1,m,n)*ts
     &         + 0.5*smodel%rlewis_coeffs(4,m,n)*t1
     &         + 1.0d0/3.0d0*smodel%rlewis_coeffs(5,m,n)*t2
     &         + 0.25*smodel%rlewis_coeffs(6,m,n)*t3
     &         + 0.2*smodel%rlewis_coeffs(7,m,n)*t4
     &         + smodel%rlewis_coeffs(9,m,n)*ti
            s  = smodel%rlewis_coeffs(3,m,n)*tl
     &         - smodel%rlewis_coeffs(2,m,n)*ti
     &         - 0.5*smodel%rlewis_coeffs(1,m,n)*ts
     &         + smodel%rlewis_coeffs(4,m,n)*t1
     &         + 0.5*smodel%rlewis_coeffs(5,m,n)*t2
     &         + 1.0d0/3.0d0*smodel%rlewis_coeffs(6,m,n)*t3
     &         + 0.25*smodel%rlewis_coeffs(7,m,n)*t4
     &         + smodel%rlewis_coeffs(10,m,n)
            cp = smodel%rlewis_coeffs(1,m,n)*ts
     &         + smodel%rlewis_coeffs(2,m,n)*ti
     &         + smodel%rlewis_coeffs(3,m,n)
     &         + smodel%rlewis_coeffs(4,m,n)*t1
     &         + smodel%rlewis_coeffs(5,m,n)*t2
     &         + smodel%rlewis_coeffs(6,m,n)*t3
     &         + smodel%rlewis_coeffs(7,m,n)*t4

            smodel%rlewis_coeffs(3,m+1,n)  = cp
            smodel%rlewis_coeffs(9,m+1,n)  = (h-cp)*t1
            smodel%rlewis_coeffs(10,m+1,n) = s - cp*tl

         endif
      enddo


c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine read_lewis_bulk_file
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ------------Function finite_rate_surf_sanity_check()-----------------
c ---------------------------------------------------------------------
      Integer Function finite_rate_surf_sanity_check(smodel)
c
c *********************************************************************
c ******    do sanity checks on the finite rate surface          ******
c ******    and issue a non-zero (C-style) error code if there   ******
c ******    is a problem (detailed warnings listed below)        ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      April    30, 2009                        ******
c ****** last updated:  May      17, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
      type(finiteRateSurfaceModelType) :: smodel
c 
c **** declaration of internally used variables
c      i,j,n,m,k  = loop counters
c      ier        = internal error flag
c      iph        = phase index (gas=1,surface=2,bulk=3)
c      iphar      = reactant phase index array
c      iphap      = product phase index array
c      irsum      = sum of reactants in each phase (gas=1,surface=2,bulk=3)
c      ipsum      = sum of products in each phase (gas=1,surface=2,bulk=3)
c      rrsum      = reactant total mass
c      rpsum      = product total mass
      Integer  :: i,j,n,m,k,iph,ier
      Integer  :: irsum(0:3),ipsum(0:3)
      Integer  :: iphar(ISMODEL_MAX_REACT),iphap(ISMODEL_MAX_PROD)
      Real*8   :: rpsum, rrsum
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN FUNCTION BODY ********************
c 

c ---- Return Values from 1 to 100 are WARNINGS
c      (continue but warn user)
c      1 = no Themodynamic data available for surface species
c      2 = pyrolysis turned on with no bulk phase
c      3 = reactants or products specified do not match the 
c          expected types for the reaction
c      4 = reaction does not balance mass
c
c ---- Return Values greater than 100 are ERRORS
c      (run should be stopped immediately!)
c      101 = An empty surface site E(s*) cannot be the first
c            reactant or product (1.0/Mw results in a divide by zero)
c
      finite_rate_surf_sanity_check = 0

c---- WARNING #1: Make sure we have desorption rates for all surface
c                 species if calculating thermo data is required
c
c     This is a WARNING, but the following occurs:
c             - ALL backward rates for reactions involving this  
c               species are implicitly turned off
c             - Gibbs energy is set to zero explicitly (kthermo = 2)
c
      k = smodel%ngps
      do n=1,smodel%nsp
        do j=1,smodel%nsps(n)
          k=k+1
	    if((smodel%kthermo(k).eq.0) .and. 
     &       (smodel%kadsr(k) .le. 0)) then

c           ---- set warning flag ----
            finite_rate_surf_sanity_check = 
     &                            max(1,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i5,a)')
     &               'missing desorption data for species ',k,
     &               ': some backward rates set to zero'
c           ---- zero Gibbs energy ----
            smodel%kthermo(k) = 2
c           ---- search through all reactions and turn off backward rates ----
            do i=1,smodel%isrt
              do m=1,ISMODEL_MAX_REACT
                if(smodel%RR(i,m) .eq. k) smodel%isrbon(i) = 0
              enddo
              do m=1,ISMODEL_MAX_PROD  
                if(smodel%RP(i,m) .eq. k) smodel%isrbon(i) = 0
              enddo
            enddo

          endif
        enddo
      enddo

c---- WARNING #2: Pyrolysis cannot be active with no bulk phase
c
      if ((smodel%isblwon .gt. 0).and.(smodel%nbp .lt. 1)) then
c        ---- set warning flag ----
         finite_rate_surf_sanity_check = 
     &                            max(2,finite_rate_surf_sanity_check)
         write(smodel%cfrm_error_string,'(a,a)')  
     &               'pyrolysis blowing cannot',
     &               ' be active without at least one bulk phase'
c        ---- turn off pyrolysis ----
         smodel%isblwon  = 0
         smodel%nblwflag = 0
      endif

c---- WARNING #3: reactions must have proper participant types
c
      do i = 1,smodel%isrt

c        ---- count the number of species in each phase ----
         irsum = 0
         ipsum = 0
         do m=1,ISMODEL_MAX_REACT
            iph = 0
            if (smodel%RR(i,m) .gt. 0) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%isend) iph = iph + 1
            irsum(iph) = irsum(iph) + smodel%RVR(i,m)
         enddo
         do m=1,ISMODEL_MAX_PROD
            iph = 0
            if (smodel%RP(i,m) .gt. 0) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%isend) iph = iph + 1
            ipsum(iph) = ipsum(iph) + smodel%RVP(i,m)
         enddo
         
c        ---- issue warning for certain types ----
         if (smodel%istype(i) .eq. 2) then      ! E-R
            if ((irsum(1) .lt. 1).or.(irsum(2) .lt. 1) .or.
     &          (ipsum(1) .lt. 1).or.(ipsum(2) .lt. 1)     ) then
               finite_rate_surf_sanity_check = 
     &                            max(3,finite_rate_surf_sanity_check)
               write(smodel%cfrm_error_string,'(a,a,i4)')  
     &               'E-R reaction ',
     &               'not of form A+[B]<-->AB+[] for reaction #',i
            endif
         elseif (smodel%istype(i) .eq. 3) then  ! L-H
            if ((irsum(1) .lt. 0).or.(irsum(2) .lt. 2) .or.
     &          (ipsum(1) .lt. 1).or.(ipsum(2) .lt. 2)     ) then
               finite_rate_surf_sanity_check = 
     &                            max(3,finite_rate_surf_sanity_check)
               write(smodel%cfrm_error_string,'(a,a,i4)')  
     &               'L-H reaction ',
     &               'not of form [A]+[B]<-->AB+2[] for reaction #',i
            endif
         elseif (smodel%istype(i) .eq. 4) then  ! sublimation
            if ((irsum(1) .lt. 0).or.(irsum(2) .lt. 1) .or.
     &          (ipsum(1) .lt. 1).or.(ipsum(2) .lt. 1) .or.
     &          ((irsum(3)+ipsum(3)) .lt. 1) ) then
               finite_rate_surf_sanity_check = 
     &                            max(3,finite_rate_surf_sanity_check)
               write(smodel%cfrm_error_string,'(a,a,i4)')  
     &               'sublimation reaction ',
     &               'not of form [A]+<B><-->AB+[] for reaction #',i
            endif
         else    !adsorption: type 1 or 5
            if ((irsum(1) .lt. 1).or.(irsum(2) .lt. 1) .or.
     &          (ipsum(1) .lt. 0).or.(ipsum(2) .lt. 1)     ) then
               finite_rate_surf_sanity_check = 
     &                            max(3,finite_rate_surf_sanity_check)
               write(smodel%cfrm_error_string,'(a,a,i4)')  
     &               'Adsorption reaction ',
     &               'not of form A+[]<-->[A] for reaction #',i
            endif
         endif
      
      enddo

c---- WARNING #4: mass balance for reactions
c
      do i = 1,smodel%isrt

c        ---- sum the mass of each reaction ----
         rrsum = 0.0d0
         rpsum = 0.0d0
         do m=1,ISMODEL_MAX_REACT
            if (smodel%RR(i,m) .gt. 0) rrsum = rrsum + 
     &          dble(smodel%RVR(i,m))*smodel%wmass(smodel%RR(i,m))
         enddo
         do m=1,ISMODEL_MAX_PROD  
            if (smodel%RP(i,m) .gt. 0) rpsum = rpsum + 
     &          dble(smodel%RVP(i,m))*smodel%wmass(smodel%RP(i,m))
         enddo
         if (rrsum .lt. 1.0d-30) rrsum = 1.0d-30
         if (rpsum .lt. 1.0d-30) rpsum = 1.0d-30
         
         if ( abs((rrsum-rpsum)/(rrsum+rpsum)) .ge. 1.0d-2 ) then
            finite_rate_surf_sanity_check = 
     &                            max(4,finite_rate_surf_sanity_check)
               write(smodel%cfrm_error_string,'(a,i4)')  
     &               'mass does not balance for reaction #',i
         endif

      enddo

c---- WARNING #5: reactions can have only one gas or bulk
c
      do i = 1,smodel%isrt

c        ---- count the number of species in each phase ----
         irsum = 0
         ipsum = 0
         do m=1,ISMODEL_MAX_REACT
            iph = 0
            if (smodel%RR(i,m) .gt. 0) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%isend) iph = iph + 1
            irsum(iph) = irsum(iph) + 1
         enddo
         do m=1,ISMODEL_MAX_PROD
            iph = 0
            if (smodel%RP(i,m) .gt. 0) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%isend) iph = iph + 1
            ipsum(iph) = ipsum(iph) + 1
         enddo
         
c        ---- issue warning for more than one gas or bulk ----
         if ((irsum(1) .gt. 1).or.(irsum(3) .gt. 1) .or.
     &       (ipsum(1) .gt. 1).or.(ipsum(3) .gt. 1)     ) then
            finite_rate_surf_sanity_check = 
     &                         max(5,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a)')  
     &            'reaction #',i,' is not an elementary reaction'
         endif
      
      enddo

c---- ERROR #101: Make sure we do not have an empty surface site
c                 as the first reactant or product
c
      do i = 1,smodel%isrt

c        ---- identify gas phase reactant participant ----      
         k = smodel%RR(i,1)
c        ---- only Arrhenius or Sublimation is allowed if wm(k) is 0.0
         if ( (smodel%istype(i) .gt. 0) .and. 
     &        (smodel%istype(i) .lt. 4) .and. 
     &        (smodel%wmass(k) .eq. 0.0) )  then

            finite_rate_surf_sanity_check = 
     &                            max(101,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,a)') 
     &               'cannot set empty ',
     &               'surface site as first reactant for types 1,2,3'
            return

         endif
c        ---- identify gas phase product participant ----      
         k = smodel%RP(i,1)
c        ---- cannot be Sublimation if wm(k) is 0.0
         if ( (smodel%istype(i) .eq. 4) .and. 
     &        (smodel%wmass(k) .eq. 0.0) )  then

            finite_rate_surf_sanity_check = 
     &                            max(101,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,a)') 
     &               'cannot set empty surface site ',
     &               'as first product for sublimation'
            return

         endif

      enddo

c---- ERROR #102: Insufficient array size for Lewis data
c
      do n = 1,smodel%nspecies
      
        if (smodel%rlewis_nrange(n) .gt. NLEWIS_FITS) then
        
           finite_rate_surf_sanity_check = 
     &                            max(102,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,a)') 
     &               'too many Lewis ranges; increase the ',
     &               'size of the NLEWIS_FITS parameter'
            return
        
        endif
      
      enddo

c---- ERROR #103: Make sure participants are in order: gas/surface/bulk
c
      do i = 1,smodel%isrt

c        ---- get phases of all participants ----
         do m=1,ISMODEL_MAX_REACT
            iph = 0
            if (smodel%RR(i,m) .gt. 0) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RR(i,m) .gt. smodel%isend) iph = iph + 1
            if (iph .eq. 0) iph = 4
            iphar(m) = iph
         enddo
         do m=1,ISMODEL_MAX_PROD  
            iph = 0
            if (smodel%RP(i,m) .gt. 0) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%igend) iph = iph + 1
            if (smodel%RP(i,m) .gt. smodel%isend) iph = iph + 1
            if (iph .eq. 0) iph = 4
            iphap(m) = iph
         enddo
c        ---- participant phases must increase
         if ( (iphar(2) .lt. iphar(1)) .or.
     &        (iphar(3) .lt. iphar(2)) .or.
     &        (iphap(2) .lt. iphap(1)) .or.
     &        (iphap(3) .lt. iphap(2)) )  then

            finite_rate_surf_sanity_check = 
     &                            max(103,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a,a)') 
     &             'reaction #',i,' participants must be ordered: ',
     &             'gas, surface, bulk, then zero entries'
            return

         endif

      enddo

c---- ERROR #104: Look for negative stoichiometric coefficients
c
      do i = 1,smodel%isrt
      
         ier = 0
         do m=1,ISMODEL_MAX_REACT
            if (smodel%RVR(i,m) .lt. 0) ier = 1
         enddo
         do m=1,ISMODEL_MAX_PROD  
            iph = 0
            if (smodel%RVP(i,m) .lt. 0) ier = 1
         enddo
c        ---- cannot be < 0
         if (ier .eq. 1)  then

            finite_rate_surf_sanity_check = 
     &                            max(104,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a,a)') 
     &             'reaction #',i,' cannot have negative ',
     &             'stoichiometric coefficients'
            return

         endif

      enddo

c---- ERROR #105: Look for bad participants
c
      do i = 1,smodel%isrt
      
         ier = 0
         do m=1,ISMODEL_MAX_REACT
            if ((smodel%RR(i,m) .lt. 0) .or. 
     &          (smodel%RR(i,m) .gt. smodel%nspecies)) ier = 1
         enddo
         do m=1,ISMODEL_MAX_PROD  
            iph = 0
            if ((smodel%RP(i,m) .lt. 0) .or. 
     &          (smodel%RP(i,m) .gt. smodel%nspecies)) ier = 1
         enddo
c        ---- cannot be < 0
         if (ier .eq. 1)  then

            finite_rate_surf_sanity_check = 
     &                            max(105,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a,a)') 
     &             'reaction #',i,' lists a bad ',
     &             'participant'
            return

         endif

      enddo

c---- ERROR #106: bad reaction types
c
      do i = 1,smodel%isrt

         if ((smodel%istype(i) .lt. 0) .or.
     &       (smodel%istype(i) .gt. 5))  then

            finite_rate_surf_sanity_check = 
     &                            max(106,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a)') 
     &             'reaction #',i,' is not a valid type'
            return

         endif

      enddo

c---- ERROR #107: bad molecular weight
c
      do k = 1,smodel%nspecies

         if (smodel%wmass(k) .lt. 0.0d0) then

            finite_rate_surf_sanity_check = 
     &                            max(107,finite_rate_surf_sanity_check)
            write(smodel%cfrm_error_string,'(a,i4,a)') 
     &          'species #',k,' cannot have negative molecular weight'
            return

         endif

      enddo

c 
c **************** END FUNCTION BODY *********************
c
      End Function finite_rate_surf_sanity_check
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ------------Subroutine finite_rate_surf_gibbs_energy()---------------
c ---------------------------------------------------------------------
      Subroutine finite_rate_surf_gibbs_energy(smodel, g, dgg, dgw, 
     &                                         tg, tw, frate, dfg, dfw)
c
c *********************************************************************
c ******    compute gibbs energies and derivatives for           ******
c ******    all species (gas, surface, and bulk)                 ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      May       4, 2009                        ******
c ****** last updated:  June     16, 2010                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c
c *** parameters
      real*8 :: PI, R, AV
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23)
c 
c *** declarations of argument types passed to the subroutine
c     smodel   = surface model with coefficient data
c     g        = gibbs array
c     dgg      = derivatives of gibbs energies w.r.t. gas tmpr, tg
c     dgw      = derivatives of gibbs energies w.r.t. wall tmpr, tw
c     tg, tw   = gas and solid state temepratures
c     frate    = array ofr forward rates for all reactions
c     dfg, dfw = derivatives of forward rates for all reactions
      type(finiteRateSurfaceModelType)  :: smodel
      real*8,dimension(smodel%nspecies) :: g, dgg, dgw
      real*8,dimension(smodel%isrt)     :: frate, dfg, dfw
      real*8  :: tg, tw
      
c 
c **** declaration of internally used variables
c     m                 = Lewis range index
c     h                 = Lewis enthalpy
c     s                 = Lewis entropy
c     dh, ds            = derivatives of Lewis fits
c     t2,t3,t4,ti,ts,tl = temperature variables
      integer :: m, n, i, kg
      real*8  :: h, s, dh, ds, eqKa, dKag, dKaw,
     &           t1, t2, t3, t4, ti, ts, tl, thd
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 
      thd = 1.0d0/3.0d0
      
      t1 = tg
      t2 = t1*t1
      t3 = t2*t1
      t4 = t3*t1
      ti = 1.0d0/t1
      ts = ti*ti
      tl = log(t1)
      
      g   = 0.0d0
      dgw = 0.0d0
      dgg = 0.0d0

c     ---- compute all gas phase species using the tg tmpr ----
      do n = 1, smodel%ngps
c       ---- determine range for this species ----
        m = 0
        do i=1,smodel%rlewis_nrange(n)
           if ( ( t1 .ge. smodel%rlewis_trange(1,i,n) ) .and. 
     &          ( t1 .le. smodel%rlewis_trange(2,i,n) ) )m = i
        enddo

        h  = smodel%rlewis_coeffs(3,m,n)
     &       + smodel%rlewis_coeffs(2,m,n)*ti*tl
     &       - smodel%rlewis_coeffs(1,m,n)*ts
     &       + 0.5*smodel%rlewis_coeffs(4,m,n)*t1
     &       + thd*smodel%rlewis_coeffs(5,m,n)*t2
     &       + 0.25*smodel%rlewis_coeffs(6,m,n)*t3
     &       + 0.2*smodel%rlewis_coeffs(7,m,n)*t4
     &       + smodel%rlewis_coeffs(9,m,n)*ti

        s  = smodel%rlewis_coeffs(3,m,n)*tl
     &       - smodel%rlewis_coeffs(2,m,n)*ti
     &       - 0.5*smodel%rlewis_coeffs(1,m,n)*ts
     &       + smodel%rlewis_coeffs(4,m,n)*t1
     &       + 0.5*smodel%rlewis_coeffs(5,m,n)*t2
     &       + thd*smodel%rlewis_coeffs(6,m,n)*t3
     &       + 0.25*smodel%rlewis_coeffs(7,m,n)*t4
     &       + smodel%rlewis_coeffs(10,m,n)

        dh = 2.0*smodel%rlewis_coeffs(1,m,n)*ts*ti
     &       + smodel%rlewis_coeffs(2,m,n)*ts*(1.0-tl)
     &       + 0.5*smodel%rlewis_coeffs(4,m,n)
     &       + 2.0*thd*smodel%rlewis_coeffs(5,m,n)*t1
     &       + 0.75*smodel%rlewis_coeffs(6,m,n)*t2
     &       + 0.8*smodel%rlewis_coeffs(7,m,n)*t3
     &       - smodel%rlewis_coeffs(9,m,n)*ts

        ds = smodel%rlewis_coeffs(1,m,n)*ts*ti
     &       + smodel%rlewis_coeffs(2,m,n)*ts
     &       + smodel%rlewis_coeffs(3,m,n)*ti
     &       + smodel%rlewis_coeffs(4,m,n)
     &       + smodel%rlewis_coeffs(5,m,n)*t1
     &       + smodel%rlewis_coeffs(6,m,n)*t2
     &       + smodel%rlewis_coeffs(7,m,n)*t3

c       FIX:  Gibbs energies are made dimensional because we have to be
c             careful with which temperature we are using!
        g(n)   = (h  - s)*R*t1
        dgg(n) = (dh - ds)*R*t1 + (h - s)*R

      enddo


c     ---- compute all surface and bulk phase species using the tw tmpr ----
      t1 = tw
      t2 = t1*t1
      t3 = t2*t1
      t4 = t3*t1
      ti = 1.0d0/t1
      ts = ti*ti
      tl = log(t1)

      do n = (smodel%ngps+1), smodel%nspecies

c       ---- we have Gibbs energy relationships for this species ----
        if ( smodel%kthermo(n) .eq. 1) then

c         ---- determine range for this species ----
          m = 0
          do i=1,smodel%rlewis_nrange(n)
             if ( ( t1 .ge. smodel%rlewis_trange(1,i,n) ) .and. 
     &            ( t1 .le. smodel%rlewis_trange(2,i,n) ) ) m = i
          enddo

          h  = smodel%rlewis_coeffs(3,m,n)
     &       + smodel%rlewis_coeffs(2,m,n)*ti*tl
     &       - smodel%rlewis_coeffs(1,m,n)*ts
     &       + 0.5*smodel%rlewis_coeffs(4,m,n)*t1
     &       + thd*smodel%rlewis_coeffs(5,m,n)*t2
     &       + 0.25*smodel%rlewis_coeffs(6,m,n)*t3
     &       + 0.2*smodel%rlewis_coeffs(7,m,n)*t4
     &       + smodel%rlewis_coeffs(9,m,n)*ti

          s  = smodel%rlewis_coeffs(3,m,n)*tl
     &       - smodel%rlewis_coeffs(2,m,n)*ti
     &       - 0.5*smodel%rlewis_coeffs(1,m,n)*ts
     &       + smodel%rlewis_coeffs(4,m,n)*t1
     &       + 0.5*smodel%rlewis_coeffs(5,m,n)*t2
     &       + thd*smodel%rlewis_coeffs(6,m,n)*t3
     &       + 0.25*smodel%rlewis_coeffs(7,m,n)*t4
     &       + smodel%rlewis_coeffs(10,m,n)

          dh = 2.0*smodel%rlewis_coeffs(1,m,n)*ts*ti
     &       + smodel%rlewis_coeffs(2,m,n)*ts*(1.0-tl)
     &       + 0.5*smodel%rlewis_coeffs(4,m,n)
     &       + 2.0*thd*smodel%rlewis_coeffs(5,m,n)*t1
     &       + 0.75*smodel%rlewis_coeffs(6,m,n)*t2
     &       + 0.8*smodel%rlewis_coeffs(7,m,n)*t3
     &       - smodel%rlewis_coeffs(9,m,n)*ts

          ds = smodel%rlewis_coeffs(1,m,n)*ts*ti
     &       + smodel%rlewis_coeffs(2,m,n)*ts
     &       + smodel%rlewis_coeffs(3,m,n)*ti
     &       + smodel%rlewis_coeffs(4,m,n)
     &       + smodel%rlewis_coeffs(5,m,n)*t1
     &       + smodel%rlewis_coeffs(6,m,n)*t2
     &       + smodel%rlewis_coeffs(7,m,n)*t3

c       FIX:  Gibbs energies are made dimensional because we have to be
c             careful with which temperature we are using!
          g(n)   = (h  - s)*R*t1
          dgw(n) = (dh - ds)*R*t1 + (h - s)*R

c       ---- trap for surface species for which there was no desorption data found ----
c            backward rates are shut-off, so the Gibbs data is unimportant
        elseif ( smodel%kthermo(n) .eq. 2 ) then

          g(n)   = 0.0d0
          dgg(n) = 0.0d0
	    dgw(n) = 0.0d0

c       ---- we need Gibbs energies extrapolated from rate data ----
        else

          i=smodel%kadsr(n)
          kg=smodel%RR(i,1)
          call adeq(tw, tg, smodel%wmass(kg), smodel%iadtype(i), 
     &              smodel%iadform(i), smodel%Cfad(i), smodel%etaad(i),
     &              smodel%vad(i),smodel%Edes(i),smodel%Ea(i),
     &              frate(i), dfg(i), dfw(i),
     &              eqKa, dKag, dKaw )

          g(n)   = g(kg) - R*tw*log(eqKa)
          dgg(n) = dgg(kg) - R*tw/eqKa*dKag
          dgw(n) = dgw(kg) - R*(log(eqKa)+tw/eqKa*dKaw)

        endif

      enddo

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine finite_rate_surf_gibbs_energy
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ---------------Subroutine finite_rate_surf_eqcon()-------------------
c ---------------------------------------------------------------------
      Subroutine finite_rate_surf_eqcon(smodel, gibbs, dgg, dgw, tw,
     &                                        EqC, dEqCg, dEqCw)
c
c *********************************************************************
c ******    compute equilibrium constants for all reactions      ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      May      14, 2009                        ******
c ****** last updated:  July      3, 2009                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c
c *** parameters
      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
c 
c *** declarations of argument types passed to the subroutine
c     smodel       = surface model with coefficient data
c     gibbs        = gibbs array
c     dgg          = derivatives of gibbs energies w.r.t. gas tmpr, tg
c     dgw          = derivatives of gibbs energies w.r.t. wall tmpr, tw
c     tw           = gas and solid state temepratures
c     EqC          = array of equilibrium constants for all reactions
c     dEqCg, dEqCw = derivatives of equilibrium constants for all reactions
      type(finiteRateSurfaceModelType)  :: smodel
      real*8,dimension(smodel%nspecies) :: gibbs, dgg, dgw
      real*8,dimension(smodel%isrt)     :: EqC, dEqCg, dEqCw
      real*8  :: tw
      
c 
c **** declaration of internally used variables
c     ireact, ispec, iphase = loop index counters
c     tmp                   = calculation variable
c     eqa                   = activity based equilibrium constant
c     deqa_dtg, deqa_dtw    = derivatives of eqa
      integer :: ireact, ispec, iphase
      real*8  :: tmp
      real*8,dimension(smodel%isrt) :: eqa, deqa_dtg, deqa_dtw
c 
c **** declaration of functions used
c      NONE


c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- initialize
      EqC = 0.0d0
      dEqCg = 0.0d0
      dEqCw = 0.0d0

      eqa      = 0.0d0
      deqa_dtg = 0.0d0
      deqa_dtw = 0.0d0
      
c---- calculate activity based equilibrium constants & derivatives
      do ireact = 1,smodel%isrt
         do ispec = 1,smodel%nspecies
            tmp = smodel%v(ireact,ispec)/(R*tw)
            eqa(ireact)      = eqa(ireact) - tmp*gibbs(ispec)
            deqa_dtg(ireact) = deqa_dtg(ireact) - tmp*dgg(ispec)
            deqa_dtw(ireact) = deqa_dtw(ireact) - tmp*
     &                    (dgw(ispec) - gibbs(ispec)/tw)
         enddo
         eqa(ireact) = exp(eqa(ireact))
         deqa_dtg(ireact) = eqa(ireact)*deqa_dtg(ireact)
         deqa_dtw(ireact) = eqa(ireact)*deqa_dtw(ireact)
      enddo


c---- calculate concentration based equilibrium constants & derivatives
      do ireact = 1,smodel%isrt
         tmp = (PREF/(R*tw))**(smodel%sumv(ireact,1))
         do iphase = 2,smodel%nsp+1
            tmp = tmp*smodel%sden(iphase-1)
     &                                  **(smodel%sumv(ireact,iphase))
         enddo
         
         eqc(ireact)   = eqa(ireact)*tmp
         dEqCg(ireact) = tmp*deqa_dtg(ireact)
         dEqCw(ireact) = tmp*(deqa_dtw(ireact) - 
     &                   smodel%sumv(ireact,1) *
     &                   eqa(ireact)/tw )
      enddo
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine finite_rate_surf_eqcon
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      Subroutine srate(Tw,Tg,wm,stot,p,itype,A,B,E,rate,dkdTg,dkdTw)
c
c     This subroutine returns the numerical value of the surface reaction
c     rate constant for the following types of reactions:
c
c     itype=0   Arrhenius 		
c     itype=1   Adsoprtion
c     itype=2   Eley-Rideal
c     itype=3   Langmuir-Hinschelwood
c     itype=4   Sublimation based on vapor pressure
c     itype=5   Arrhenius adsorption
c
c     written by: Jochen Marschall
c     updated by: Matthew MacLean (maclean@cubrc.org)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Implicit none

      real*8 :: PI, R, AV
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23)
c 
c *** declarations of argument types passed to the subroutine
c     Tw      = solid state surface temperature
c     Tg      = near-wall gas phase temerature
c     wm      = molecular weight
c     stot    = total surface density
c     p       = site density exponent
c     itype   = type of surface reaction
c     A, B, E = reaction specific parameters
c     rate    = forward rate to return
c     dkdTg   = derivative w.r.t. gas temperature
c     dkdTw   = derivative w.r.t. wall temperature
      integer :: itype
      real*8  :: Tw,Tg,wm,stot,p,A,B,E,rate,dkdTg,dkdTw
c 
c **** declaration of internally used variables
      real*8 :: v,S0c,S0,gERc,gER,v2D,clh,a0vp0

      dkdTg = 0.0d0
      dkdTw = 0.0d0

c-----Arrhenius Reaction
      if((itype.eq.0).or.(itype.eq.5)) then
        rate  = A*(Tw**B)*exp(-E/(Tw*R))
        dkdTw = rate * (B + E/(Tw*R))/Tw

c-----Adsorption Reaction: A=sticking coefficient, E=Adsorption Energy barrier
      else if (itype.eq.1) then
        v=sqrt(8.0*R*Tg/(PI*wm))
        S0c=A*Tw**B
        S0=min(1.0d0,S0c)
        rate=(v/(4.0*(stot**p)))*S0*exp(-E/(Tw*R))
        if (S0c .gt. 1.0d0) then
           dkdTw = rate * E/(Tw*R*Tw)
        else
           dkdTw = rate * (B+E/(Tw*R))/Tw
        endif
        dkdTg = rate * 0.5d0/Tg

c-----Eley-Rideal Reaction: A=reaction efficiency, E=Reaction Energy barrier
      else if (itype.eq.2) then
        v=sqrt(8.0*R*Tg/(PI*wm))
        gERc=A*Tw**B
        gER=min(1.0d0,gERc)
        rate=(v/(4.0*(stot**p)))*gER*exp(-E/(Tw*R))
        if (gERc .gt. 1.0d0) then
           dkdTw = rate * E/(Tw*Tw*R)
        else
           dkdTw = rate * (B+E/(Tw*R))/Tw
        endif
        dkdTg = rate * 0.5d0/Tg

c-----Langmuir-Hinschelwood Reaction: A=reaction efficiency, E=Reaction barrier
      else if (itype.eq.3) then
        v2D=sqrt(PI*R*Tw/(2.0*wm))
        clh = A*Tw**B
        rate=v2D*sqrt(AV*stot**(3.0-2.0*p))*clh*exp(-E/(Tw*R))
        dkdTw = rate * (0.5d0+B+E/(Tw*R))/Tw

c-----Sublimation Reaction: based on vapor pressure
      else if (itype.eq.4) then
        v=sqrt(8.0*R*Tg/(PI*wm))
        a0vp0=(A*Tw**B)/(R*Tg)
        rate=(v/(4.0*stot**p))*a0vp0*exp(-E/(Tw*R))
        dkdTw = rate * (B + E/(Tw*R))/Tw
        dkdTg = -0.5*rate/Tg

      endif

      return
      end Subroutine srate

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      Subroutine adeq(Tw,Tg,wm,itype,iform,A,B,v,Ed,Ea,frate,dfdTg,
     &                dfdTw,eqKa,dEQdTg,dEQdTw)
c
c      This subroutine returns the activity-based eqilibrium constant
c      for the adsorbed species at temperature Tg/Tw using the specified
c      backwards reaction or equilibrium constant parameters 
c      and the temperatures
c
c      Type 1  Desorption
c  	  Form 0:  Arrhenius
c	  Form 1:  constant frequency
c       Form 2:  Simple TST frequency
c       Form 3:  Complex TST frequency
c
c      Type 2  Equilibrium constant
c	  Form 0:  Arrhenius
c       Form 1:  Immobile simple TST
c       Form 2:  Immobile complex TST
c       Form 3:  Mobile simple TST
c       Form 4:  Mobile complex TST
c
c     written by: Jochen Marschall
c     updated by: Matthew MacLean (maclean@cubrc.org)
c         
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      
      real*8 :: PI, R, AV, Pref, H
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23)
      parameter(Pref=1.0e5,H=6.62618e-34)
c 
c *** declarations of argument types passed to the subroutine
c     Tw          = solid state surface temperature
c     Tg          = near-wall gas phase temerature
c     wm          = molecular weight
c     stot        = total surface density
c     itype       = type of surface reaction
c     A,B,v,Ed,Ea = reaction specific parameters
c     itype       = type of reaction (backward rate or equilibrium)
c     iform       = form of equilibrium or backward rate relation
c     frate       = existing forward rate
c     dfdTg,dfdTw = forward rate derivatives
c     eqKa        = value of equilibrium constant to return
c     dEQdTg      = derivative w.r.t. gas temperature
c     dEQdTw      = derivative w.r.t. wall temperature      
      Integer :: itype, iform
      real*8  :: Tw,Tg,wm,A,B,v,Ed,Ea,frate,eqKa,dEQdTg,dEQdTw,
     &           dfdTg,dfdTw

c 
c **** declaration of internally used variables
      real*8 :: brate, vfac, eqKc, qtr, qvib, qtr1, qtr2, dkbdTw,
     &          dvf, eterm

      dEQdTg = 0.0d0
      dEQdTw = 0.0d0
      dkbdTw = 0.0d0

c--------------------------------------------------`
c-----Desorption reaction
      if (itype.eq.1) then

c------ Form 0:  Arrhenius
        if (iform.eq.0) then
          brate=A*(Tw**B)*exp(-Ed/(Tw*R))
          dkbdTw = brate * (B + Ed/(R*Tw))/Tw

c------ Form 1: Constant attempt frequency
        else if (iform.eq.1) then
          brate=A*(Tw**B)*v*exp(-Ed/(Tw*R))
          dkbdTw = brate * (B + Ed/(R*Tw))/Tw

c------ Form 2: Simple TST frequency
        else if (iform.eq.2) then
          brate=A*(Tw**B)*((R*Tw)/(AV*H))*exp(-Ed/(Tw*R))
          dkbdTw = brate * (B + 1.0d0 + Ed/(R*Tw))/Tw

c------ Form 3: complex TST frequency
        else if (iform.eq.3) then
          eterm = exp((R*Tw)/(Av*H*v))
          vfac=sqrt(eterm)/(1.0-eterm)
          brate=A*(Tw**B)*((R*Tw)/(AV*H))/vfac*exp(-Ed/(Tw*R))
          dvf = 0.5d0*R/(Av*H*v)*(1.0+eterm)/(1.0-eterm)
          dkbdTw = brate * ((B + 1.0d0 + Ed/(R*Tw))/Tw - dvf)
        endif

        eqKc=frate/brate
        dEQdTw = dfdTw/brate - frate*dkbdTw/(brate*brate)
        dEQdTg = dfdTg/brate

c-------------------------------------------------------
c-----Equilibrium constant
      else if (itype.eq.2) then

c------ Form 0:  Arrehenius   
        if (iform.eq.0) then
          eqKc=A*(Tw**B)*exp(Ed/(Tw*R))
          dEQdTw = eqKc * (B - Ed/(Tw*R))/Tw

c------ Form 1:  Immobile simple TST
        else if (iform.eq.1) then
          qtr=( (2.0*PI*wm*R*Tg)/(AV*H*AV*H) )**(-1.5)
          eqKc=A*(Tw**B)*qtr*exp(-Ea/(Tw*R))*exp(Ed/(Tw*R))
          dEQdTw = eqKc * (B + (Ea-Ed)/(Tw*R))/Tw
          dEQdTg = -1.5d0*eqKc/Tg

c------ Form 2:  Immobile complex TST
        else if (iform.eq.2) then
          eterm = exp((-H*AV*v)/(R*Tw))
          qtr=( (2.0*PI*wm*R*Tg)/(AV*H*AV*H) )**(-1.5)
          qvib=( sqrt(eterm)/(1.0-eterm) )**3.0
          eqKc=A*(Tw**B)*qtr*qvib*exp(-Ea/(Tw*R))*exp(Ed/(Tw*R))
          dEQdTw = eqKc * (B + (Ea-Ed)/(Tw*R) + 1.5d0*
     &               (H*AV*v/(R*Tw))*(1.0d0+eterm)/(1.0d0-eterm))/Tw
          dEQdTg = -1.5d0*eqKc/Tg

c------ Form 3:  Mobile simple TST
        else if (iform.eq.3) then
          qtr1=( (2.0*PI*wm*R*Tg)/(AV*H*AV*H) )**(-1.5)
          qtr2=( (2.0*PI*wm*R*Tw)/(AV*H*AV*H) )
          eqKc=A*(Tw**B)*qtr1*qtr2*exp(-Ea/(Tw*R))*exp(Ed/(Tw*R))
          dEQdTw = eqKc * (B + 1.0d0 + (Ea-Ed)/(Tw*R))/Tw
          dEQdTg = -1.5d0*eqKc/Tg

c------ Form 4:  Mobile complex TST
        else if (iform.eq.4) then
          qtr1=( (2.0*PI*wm*R*Tg)/(AV*H*AV*H) )**(-1.5)
          qtr2=( (2.0*PI*wm*R*Tw)/(AV*H*AV*H) )
          eterm = exp((-H*AV*v)/(R*Tw))
          qvib= ( sqrt(eterm)/(1.0-eterm) ) 
          eqKc=A*(Tw**B)*qtr1*qtr2*qvib*exp(-Ea/(Tw*R))*exp(Ed/(Tw*R))
          dEQdTw = eqKc * (B + 1.0d0 + (Ea-Ed)/(Tw*R) + 0.5d0*
     &               (H*AV*v/(R*Tw))*(1.0d0+eterm)/(1.0d0-eterm))/Tw
          dEQdTg = -1.5d0*eqKc/Tg
        
        endif

      endif

c-----Calculate activity-based equilibrium constant
c      eqKa=eqKc*(R*Tg/Pref)
c      dEQdTw = dEQdTw*(R*Tg/Pref)
c      dEQdTg = dEQdTg*(R*Tg/Pref) + eqKc*(R/Pref)
      eqKa=eqKc*Pref/(R*Tg)
      dEQdTw = dEQdTw*Pref/(R*Tg)
      dEQdTg = dEQdTg*Pref/(R*Tg) - eqKc*Pref/(R*Tg*Tg)

      return
      end Subroutine adeq
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
