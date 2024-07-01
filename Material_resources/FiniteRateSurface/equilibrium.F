c --------------------------------------------------------------------- 
c --------------Subroutine srmodel_compute_diagnostics()---------------
c ---------------------------------------------------------------------
      Subroutine srmodel_compute_diagnostics(cell)
c
c *********************************************************************
c ******    output diagnostic rates for all reactions            ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      July     12, 2009                        ******
c ****** last updated:  April    28, 2011                        ******
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
      type(finiteRateSurfaceCell)                :: cell
c 
c **** declaration of internally used variables
c     ier           = error flag for subroutine
c     ireact        = reaction loop counter
c     ispec         = species loop counter
c     k             = misc. integer
c     pf,pb,eqact   = actual forward & backward rate source terms
c     eqact         = actual EqC computed from current concentrations
c     smodel        = copy of pointer to the surface model
c     frates,brates = kf and kb values all reactions
c     eqc           = computed equilibrium constants
c     gibbs         = list of gibbs energies
c     dfg, dfw      = forward rate derivatives (ignored)
c     dgg, dgw      = Gibbs energy derivatives (ignored)
      Integer :: ireact, ispec, k
      real*8  :: pf, pb, eqact
      type(finiteRateSurfaceModelType),pointer :: smodel
      real*8,dimension(cell%frs_record%isrt) :: frates, brates, eqc, 
     &                                          dfg, dfw
      real*8,dimension(cell%frs_record%nspecies) :: gibbs, dgw, dgg

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- copy model pointer (easy access)
      smodel => cell%frs_record
      

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
     &              dfg(ireact), dfw(ireact))
   
      enddo


c---- obtain equilibrium constants
      call finite_rate_surf_gibbs_energy(smodel, gibbs, dgg, dgw,
     &                                   cell%gas_tmpr, cell%wall_tmpr,
     &                                   frates, dfg, dfw)


      call finite_rate_surf_eqcon(smodel, gibbs, dgg, dgw,
     &                   cell%wall_tmpr, eqc, dfg, dfw)
      
      
c---- compute backward rates
      do ireact = 1,smodel%isrt
         brates(ireact) = frates(ireact)/eqc(ireact)
      enddo


c---- output diagnostics      
      write(6,*)
      write(6,600)
      write(6,601) 'Reaction Constants'
      write(6,602)
      cell%srhs(:) = 0.0d0
      do ireact = 1,smodel%isrt

c        ----- sum forward and backward concentration product term -----
         pf    = smodel%isrfon(ireact)*frates(ireact)*
     &           smodel%sfraction(ireact)
         pb    = smodel%isrbon(ireact)*brates(ireact)*
     &           smodel%sfraction(ireact)
         eqact = 1.0d0
         do ispec = 1,smodel%nspecies
            pf = pf*cell%concentration(ispec)**smodel%vr(ireact,ispec)
            pb = pb*cell%concentration(ispec)**smodel%vp(ireact,ispec)
            eqact = eqact * 
     &              cell%concentration(ispec)**smodel%v(ireact,ispec)
         enddo
         
         do ispec = 1,smodel%nspecies
            cell%srhs(ispec) = cell%srhs(ispec) + 
     &                         smodel%v(ireact,ispec) * (pf - pb)
         enddo

c        ----- write reaction diagnostic statement -----
         write(6,603) ireact,frates(ireact),brates(ireact),eqc(ireact),
     &                eqact,pf,pb

      enddo

c     ----- write species production diagnostic statement -----
      write(6,601) 'Species Production'
      write(6,604) 
      do ispec = 1,smodel%nspecies
         write(6,605) smodel%species(ispec),
     &                cell%concentration(ispec),cell%srhs(ispec),
     &                smodel%wmass(ispec)*cell%srhs(ispec)
      enddo

c     ----- write species production per reaction -----
      write(6,601) 'Species Production Per Reaction'
      write(6,606)
      do ireact = 1,smodel%isrt
         if ((ireact .lt. 10))   write(6,607) ireact
         if ((ireact .ge. 10).and.(ireact .lt. 100))  
     &                                             write(6,608) ireact
         if ((ireact .ge. 100).and.(ireact .lt. 1000))
     &                                             write(6,609) ireact
      enddo
      write(6,*)

      do ispec = 1,smodel%nspecies
         write(6,610) smodel%species(ispec),cell%srhs(ispec)
         do ireact = 1,smodel%isrt
            pf    = smodel%isrfon(ireact)*frates(ireact)*
     &              smodel%sfraction(ireact)
            pb    = smodel%isrbon(ireact)*brates(ireact)*
     &              smodel%sfraction(ireact)
            do k = 1,smodel%nspecies
              pf = pf*
     &             cell%concentration(k)**smodel%vr(ireact,k)
              pb = pb*
     &             cell%concentration(k)**smodel%vp(ireact,k)
            enddo

            write(6,611) smodel%v(ireact,ispec) * (pf - pb)
         enddo
         write(6,*)
      enddo
      


 600  format('---Additional Diagnostics-----------------------',
     &       '------------------------------------')
 601  format('      <<',a,'>>')
 602  format('      react      kf          kb         KeqC    ',
     &       ' KeqC-actual   wf-actual   wb-actual')
 603  format('       ',i3,2x,6es12.4)
 604  format('      species      concentration     ',
     &       'molar prod.      mass prod.')
 605  format('      ',a,3es16.4)
 606  format('      species          w-dot',$)
 607  format('             w',i1,$)
 608  format('            w',i2,$)
 609  format('           w',i3,$)
 610  format('      ',a,es12.4,$)
 611  format('   ',es12.4,$)
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine srmodel_compute_diagnostics
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ---------------Subroutine srmodel_scan_diagnostics()-----------------
c ---------------------------------------------------------------------
      Subroutine srmodel_scan_diagnostics(cell)
c
c *********************************************************************
c ******    output scanmode production branching                 ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May       2, 2011                        ******
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
c *** parameters
      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
c 
c *** declarations of argument types passed to the subroutine
c     cell = cell with current concentrations already set
      type(finiteRateSurfaceCell)                :: cell
c 
c **** declaration of internally used variables
c     ier           = error flag for subroutine
c     ireact        = reaction loop counter
c     ispec         = species loop counter
c     k             = misc. integer
c     pf,pb         = actual forward & backward rate source terms
c     srhs          = species total production
c     smodel        = copy of pointer to the surface model
c     frates,brates = kf and kb values all reactions
c     eqc           = computed equilibrium constants
c     gibbs         = list of gibbs energies
c     dfg, dfw      = forward rate derivatives (ignored)
c     dgg, dgw      = Gibbs energy derivatives (ignored)
      Integer :: ireact, ispec, k
      real*8  :: pf, pb, srhs
      type(finiteRateSurfaceModelType),pointer :: smodel
      real*8,dimension(cell%frs_record%isrt) :: frates, brates, eqc, 
     &                                          dfg, dfw
      real*8,dimension(cell%frs_record%nspecies) :: gibbs, dgw, dgg

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- copy model pointer (easy access)
      smodel => cell%frs_record
      

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
     &              dfg(ireact), dfw(ireact))
   
      enddo


c---- obtain equilibrium constants
      call finite_rate_surf_gibbs_energy(smodel, gibbs, dgg, dgw,
     &                                   cell%gas_tmpr, cell%wall_tmpr,
     &                                   frates, dfg, dfw)


      call finite_rate_surf_eqcon(smodel, gibbs, dgg, dgw,
     &                   cell%wall_tmpr, eqc, dfg, dfw)
      
      
c---- compute backward rates
      do ireact = 1,smodel%isrt
         brates(ireact) = frates(ireact)/eqc(ireact)
      enddo


c---- output diagnostics: species production per reaction
      srhs = 0.0d0
      do ispec = 1,smodel%nspecies
         do ireact = 1,smodel%isrt
            pf    = smodel%isrfon(ireact)*frates(ireact)*
     &              smodel%sfraction(ireact)
            pb    = smodel%isrbon(ireact)*brates(ireact)*
     &              smodel%sfraction(ireact)
            do k = 1,smodel%nspecies
              pf = pf*
     &             cell%concentration(k)**smodel%vr(ireact,k)
              pb = pb*
     &             cell%concentration(k)**smodel%vp(ireact,k)
            enddo
            write(6,701) smodel%v(ireact,ispec) * (pf - pb)
            srhs = srhs + smodel%v(ireact,ispec) * (pf - pb)
         enddo
      enddo

 701  format(x,es14.6,$)
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine srmodel_scan_diagnostics
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c --------------Subroutine srmodel_compute_equilibrium()---------------
c ---------------------------------------------------------------------
      Subroutine srmodel_compute_equilibrium(cell, eta, igastype)
c
c *********************************************************************
c ******    free-energy minimization wrapper                     ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      June     22, 2009                        ******
c ****** last updated:  October  12, 2010                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none

c
c *** interface declaration required 
      interface
         Subroutine generate_element_table(nspecies,snames,
     &                                     etable,nelements,ierror)
            Integer :: nspecies, nelements, ierror
            Character*10,dimension(nspecies)   :: snames
            Real*8,dimension(:,:),pointer     :: etable
         end Subroutine generate_element_table
      end interface

c
c *** parameters
      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
c 
c *** declarations of argument types passed to the subroutine
c     cell     = cell with current concentrations already set
c     eta      = array of initial values of mole-mass and conc
c     igastype = flag for gas phase treatment
      type(finiteRateSurfaceCell)                :: cell
      real*8,dimension(cell%frs_record%nspecies) :: eta
      integer                                    :: igastype
c 
c **** declaration of internally used variables
c     ier           = error flag for subroutine
c     ireact        = reaction loop counter
c     iphase        = phase associated with a surface species
c     k             = misc. integer
c     res           = residual from solver
c     smodel        = copy of pointer to the surface model
c     dfg, dfw      = forward rate derivatives (ignored)
c     dgg, dgw      = Gibbs energy derivatives (ignored)
c     frates        = forward rates
c     gibbs         = Gibbs values
c     sref          = reference site density
      Integer :: ier, ireact, k, iphase
      real*8  :: res
      type(finiteRateSurfaceModelType),pointer :: smodel
      real*8,dimension(cell%frs_record%isrt) :: frates, dfg, dfw
      real*8,dimension(cell%frs_record%nspecies) :: gibbs, dgw, dgg, 
     &                                              sref

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 

c---- copy model pointer (easy access)
      smodel => cell%frs_record
      

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
     &              dfg(ireact), dfw(ireact))
   
      enddo


c---- obtain Gibbs energies for all species either by Lewis data or backward rates
      call finite_rate_surf_gibbs_energy(smodel, gibbs, dgg, dgw,
     &                                   cell%gas_tmpr, cell%wall_tmpr,
     &                                   frates, dfg, dfw)


c---- Gibbs energies need to be non-dimensional for this to work
      do k=1,cell%num_gas
         gibbs(k) = gibbs(k) / (R*cell%gas_tmpr)
      enddo
      do k=(cell%num_gas+1),cell%num_species
         gibbs(k) = gibbs(k) / (R*cell%wall_tmpr)
      enddo


c---- get reference site densities for the surface species
      sref = 1.0d0    !make ln(sref)=0.0 by default
      do k=(cell%num_gas+1),cell%num_nonbulk
         iphase  = smodel%kphase(k)
         sref(k) = smodel%sden(iphase-1)
      enddo


c---- call the solver
      call free_energy_minimization(cell%num_species, cell%num_gas, 
     &                              cell%num_nonbulk, igastype,
     &                              smodel%species, gibbs, eta, 
     &                              cell%wall_pres, cell%gas_tmpr, 
     &                              sref, res, ier)
c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine srmodel_compute_equilibrium
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c --------------Subroutine free_energy_minimization()------------------
c ---------------------------------------------------------------------
      Subroutine free_energy_minimization(nspecies, ngas, nonbulk, 
     &                                    igastype, snames, gibbs, eta,
     &                                    pres, tmpr, sref, res, 
     &                                    ierror)
c
c *********************************************************************
c ******    perform a free-energy minimization calculation       ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      June     22, 2009                        ******
c ****** last updated:  May      17, 2011                        ******
c *********************************************************************
c
      Implicit none

c     ----- interface for element table subroutine -----
      interface
         Subroutine generate_element_table(nspecies,snames,
     &                                     etable,nelements,ierror)
            Integer :: nspecies, nelements, ierror
            Character*10,dimension(nspecies)  :: snames
            Real*8,dimension(:,:),pointer     :: etable
         end Subroutine generate_element_table
      end interface


c     ----- declaration of external arguments -----
      Integer :: nspecies, ngas, nonbulk, nvol, ierror, igastype
      Real*8  :: pres, tmpr, res
      Character*10,dimension(nspecies)   :: snames
      Real*8,dimension(nspecies)         :: gibbs, eta, sref


c     ----- declaration of internally supplied variables -----
      Real*8,dimension(:,:),pointer     :: etable
      Real*8,allocatable,dimension(:)   :: uval, rhs, 
     &                                     bi0_array
      Real*8,allocatable,dimension(:,:) :: rjacobian,ebal
      Integer,allocatable,dimension(:)  :: ibasis
      Integer,dimension(nspecies)       :: inc_species
      Real*8,dimension(nspecies)        :: gtot
      Real*8  :: tterm, rtmp, volume
      Integer :: i, j, k, index, icount, nsurf, nbulk, nelements,
     &           neqns, nfrozen, idx, ilen, iout, iin
c 
c **************** BEGIN SUBROUTINE BODY *******************
c
      nsurf = nonbulk - ngas
      nbulk = nspecies - nonbulk

      nfrozen = nbulk
      if (igastype.eq.1) nfrozen = nfrozen + ngas

      nvol = 0
      if (igastype.eq.2) nvol = 1
      volume = 1.0d0
      
      inc_species = 1
      
      gtot = 0.0d0
      
      tterm = log(8.314d0*tmpr/1.0d5)

c     ----- FIRST: Generate a table of elemental constraints -----
      call generate_element_table(nspecies, snames, 
     &                            etable, nelements, ierror)

c     ----- SECOND: catch frozen redundancies ---
      if (nfrozen .gt. 0) then
         allocate(ebal(nelements,nfrozen),stat=ierror)
         allocate(ibasis(nelements),stat=ierror)
         idx = 1
         if (igastype.eq.1) then
            do j = 1,ngas
               inc_species(j) = 0
               rtmp = max(eta(j),1.0d-20)
               gtot(idx) = gibbs(j) + log(rtmp) + tterm
               do i = 1,nelements
                  ebal(i,idx) = etable(j,i)
               enddo
               idx = idx + 1
            enddo
         endif
         do j = 1,nbulk
            inc_species(nonbulk+j) = 0
            gtot(idx) = gibbs(nonbulk+j)
            do i = 1,nelements
               ebal(i,idx) = etable(nonbulk+j,i)
            enddo
            idx = idx + 1
         enddo

c        -- use row/echelon to select species --
c         call rowechelon(ebal, nelements, nfrozen, ibasis)
         
c        -- use simplex to select species --         
         call simplex(ebal, nelements, nfrozen, gtot, ibasis)

c        -- set active species --
         do i = 1,nelements
            j = ibasis(i)
            if (j .gt. 0) then
               if (igastype.eq.1) then
                  if (j .gt. ngas) j = j - ngas + nonbulk
               else
                  j = j + nonbulk
               endif
               inc_species(j) = 1
            endif
         enddo

c        -- write active system --
         write(6,*) '     Gibbs Data of Species:'
         do j = 1,nspecies
            rtmp = max(eta(j),1.0d-20)
            if (j .le. ngas) then
               rtmp = gibbs(j) + log(rtmp) + tterm
            else
               rtmp = gibbs(j)
            endif
            write(6,'(a,a,2es14.4)') '       ',
     &                            snames(j),gibbs(j),rtmp
         enddo
         write(6,*) '     Species Considered for System:'
         write(6,'(a,$)') '       '
         do j = 1,nspecies
            ilen = len_trim(snames(j))
            if (inc_species(j) .gt. 0) 
     &                      write(6,'(a,a,$)') '  ',snames(j)(1:ilen)
         enddo
         write(6,*)

         deallocate(ebal,ibasis,stat=ierror)
      endif

c     ----- THIRD: generate arrays -----
      neqns    = nspecies + nelements + nvol
      allocate(bi0_array(nelements),stat=ierror)
      allocate(uval(neqns),rhs(neqns),
     &         rjacobian(neqns,neqns),stat=ierror)


c     ----- FOURTH: initialize field -----
      uval      = 0.0d0
      rhs       = 0.0d0
      rjacobian = 0.0d0
      bi0_array = 0.0d0
      
      do i = 1,nelements
         do j = 1,nspecies
            bi0_array(i) = bi0_array(i) + etable(j,i)*eta(j)
         enddo
      enddo
      
      uval = 0.0d0
      do k = 1,nspecies
         uval(k) = max(eta(k),1.0d-7)
      enddo

c     ----- MAIN LOOP -----
      do icount = 1,100
      
         rjacobian = 0.0d0
         rhs       = 0.0d0

c        --- Gas phase Gibbs equations ---
         do j = 1, ngas
            
            if (igastype.ne.1) rjacobian(j,j) =  1.0d0

            do i = 1,nelements
               index = nspecies + i
               rjacobian(j,index) = etable(j,i)
            enddo
            
            rhs(j) = -(gibbs(j) + log(uval(j)) + tterm)
         
         enddo

c        --- Surface phase Gibbs equations ---
         do j = ngas + 1,ngas + nsurf

            rjacobian(j,j) =  1.0d0

            do i = 1,nelements
               index = nspecies + i
               rjacobian(j,index) = etable(j,i)
            enddo
            
            rhs(j) = -( gibbs(j) + log(uval(j)) + log(sref(j)) )

         enddo

c        --- Bulk phase Gibbs equations ---
         do j = nonbulk+1,nspecies

            do i = 1,nelements
               index = nspecies + i
               rjacobian(j,index) = etable(j,i)
            enddo

            rhs(j) = -( gibbs(j) + log(eta(j)) )

         enddo

c        --- elemental constraints ---
         do i = 1,nelements
            index = nspecies + i

            rhs(index) = bi0_array(i)
            do j = 1,ngas
               rjacobian(index,j) = volume*etable(j,i)*uval(j)
               rhs(index) = rhs(index) - volume*etable(j,i)*uval(j)
            enddo
            do j = ngas+1,nspecies
               rjacobian(index,j) = etable(j,i)*uval(j)
               rhs(index) = rhs(index) - etable(j,i)*uval(j)
            enddo
            if(igastype.eq.2) then
               do j = 1,ngas
                  rjacobian(index,neqns) = rjacobian(index,neqns)
     &                                   + etable(j,i)*uval(j)
               enddo
            endif
            
         enddo

c        --- volume contraint (if constant pressure) ---
         if (igastype.eq.2) then
            rhs(neqns) = pres/(8.314d0*tmpr)
            do j = 1,ngas
               rhs(neqns) = rhs(neqns) - uval(j)
               rjacobian(neqns,j) = uval(j)
            enddo
         endif

c        --- force out non-participating species ---
         do j = 1,nspecies
            if (inc_species(j) .eq. 0) then
               do i = 1,neqns
                  rjacobian(j,i) = 0.0d0
               enddo
               rjacobian(j,j) = 1.0d0
               rhs(j) = 0.0d0
            endif
         enddo

c         do i = 1,neqns
c            do j = 1,neqns
c               write(6,200) rjacobian(i,j)
c            enddo
c            write(6,201) rhs(i)
c         enddo
c 200     format(es10.2,' ',$)
c 201     format('|| ',es10.2)

c        --- solve ---
         call lusolve(rjacobian, rhs, neqns)

c        --- limit & update ---
         do k = 1,nspecies
            if (rhs(k) .gt. 1.0d0)  rhs(k) =  1.0d0
            if (rhs(k) .lt. -1.0d0) rhs(k) = -1.0d0
         enddo
         
         do j = 1, nspecies
            uval(j) = exp(log(uval(j)) + rhs(j))
         enddo

      enddo

c     ----- LAST: Store the values for return -----
      res = 0.0d0
      do k = 1,nspecies
         if (uval(k) .gt. 1.0d-15) res = res + rhs(k)*rhs(k)
      enddo
      
      if (igastype.ne.1) then
         do k = 1,ngas
            eta(k) = uval(k)
         enddo
      endif
      do k = ngas+1,nonbulk
         eta(k) = uval(k)
      enddo

c     ----- FINALLY: deallocate memory & cleanup -----      
      if (associated(etable))   deallocate(etable)
      
      if (allocated(bi0_array)) deallocate(bi0_array)
      if (allocated(uval))      deallocate(uval)
      if (allocated(rhs))       deallocate(rhs)
      if (allocated(rjacobian)) deallocate(rjacobian)

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine free_energy_minimization
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ----------------------Subroutine simplex()---------------------------
c ---------------------------------------------------------------------
      Subroutine simplex(ebal, nelements, nspecies, gibbs, ibasis)
c
c *********************************************************************
c ******    perform simplex method on elemental constraints      ******
c ******    to obtain minimum energy choices                     ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      May      14, 2011                        ******
c ****** last updated:  May      17, 2011                        ******
c *********************************************************************
c

c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     nelements = number of nelements
c     nspecies  = number of species
c     gibbs     = list of Gibbs energies
c     ebal      = matrix of elemental constraints
c     ibasis    = list of basis vector for each row
      Integer                               :: nelements, nspecies
      real*8,dimension(nspecies)            :: gibbs
      real*8,dimension(nelements, nspecies) :: ebal
      Integer,dimension(nelements)          :: ibasis
c 
c *** declaration of internally used variables
c     rsys        = matrix of optimization system
c     nrows,ncols = size of system
c     i, j, k     = counters
c     iprow,ipcol = pivot row and column
c     tmp         = temporary storage
      real*8,allocatable,dimension(:,:) :: rsys
      integer,dimension(nelements) :: isactive
      Integer :: nrows, ncols, i, j, k, iprow, ipcol, nactive
      real*8  :: tmp
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN SUBROUTINE BODY ********************
c 	
      ibasis = 0
	
      isactive = 0
      nactive = 0
      do i = 1,nelements
         do j = 1,nspecies
            if (ebal(i,j) .gt. 0) isactive(i) = 1
         enddo
         if (isactive(i) .gt. 0) nactive = nactive + 1
      enddo
      
      nrows = nactive + 1 + nspecies + 2
      ncols = nspecies + 2*nactive + 1 + nspecies + 1
      allocate(rsys(nrows,ncols))
      rsys = 0.0d0
      
c      write(6,*) 'nelements = ',nelements,' nspecies = ',nspecies
c      write(6,*) 'nactive = ',nactive

c     ----- build system with constraints -----

c        --- constraints: elemental composition ---
      k = 0
      do i = 1,nelements
         if (isactive(i) .gt. 0) then
            k = k + 1
            do  j = 1,nspecies
               if (ebal(i,j) .gt. 0) rsys(k,j) = 1.0d0
            enddo
            rsys(k,nspecies+k)           = -1.0d0
            rsys(k,2*nspecies+nactive+k) = 1.0d0
            rsys(k,ncols)                = 1.0d0
         endif
      enddo
c        --- constraint: sum[xj] = Ne ---
      do j = 1,nspecies
         rsys(nactive+1,j) = 1.0d0
      enddo
      rsys(nactive+1,ncols-1) = 1.0d0
      rsys(nactive+1,ncols)   = dble(nactive)
c        --- constraints: xj <= 1 ---
      do j = 1,nspecies
         rsys(nactive+1+j,j) = 1.0d0
         rsys(nactive+1+j,nspecies+nactive+j) = 1.0d0
         rsys(nactive+1+j,ncols) = 1.0d0
      enddo
c        --- actual cost function: sum[gibbs(j)*xj] ---
      tmp = 0.0d0
      do j = 1,nspecies
         if (tmp .lt. gibbs(j)) tmp = gibbs(j)
c         write(6,*) 'got Gibbs  = ',gibbs(j)
      enddo
c      write(6,*) 'found MAX Gibbs  = ',tmp
      do j = 1,nspecies
         rsys(nspecies+nactive+2,j) = gibbs(j) - tmp - 1.0d0
      enddo
c        --- artificial cost function: sum[zk] ---
      do i = 1,nactive+1
         do j=1,nspecies+nactive
            rsys(nrows,j) = rsys(nrows,j) - rsys(i,j)
         enddo
         rsys(nrows,ncols) = rsys(nrows,ncols) - rsys(i,ncols)
      enddo
      
c      write(6,*) ' -------------------------------------- '
c      do i = 1,nrows
c         write(6,'(20f6.2)') rsys(i,:)
c      enddo

c     ----- perform phase I simplex -----
      do k = 1,10*ncols
         ipcol = 0
         iprow = 0
         tmp = 0.0d0

c        --- pivot column is most negative artificial cost function ---
         do j = 1,ncols-1
            if (rsys(nrows,j) .lt. tmp) then
               tmp = rsys(nrows,j)
               ipcol = j
            endif
         enddo

         if (ipcol .eq. 0) goto 201

c        --- pivot row has the smallest positive ratio ---
         tmp = 0.0d0
         do j = 1,nrows - 2
            if ( (rsys(j,ipcol) .gt. 0.0d0) .and. 
     &           (rsys(j,ncols) .ge. 0.0d0) ) then
               if ((tmp .gt. (rsys(j,ncols) / rsys(j,ipcol))) .or.
     &              (iprow .eq. 0)) then
                  tmp = rsys(j,ncols) / rsys(j,ipcol)
                  iprow = j
c                  write(6,*) 'found j=',j,' tmp=',tmp,' iprow=',iprow
               endif
            endif
         enddo
         
         if (iprow .eq. 0) goto 201

c        --- normalize pivot row ---
         tmp = rsys(iprow,ipcol)
         do j = 1,ncols
            rsys(iprow,j) = rsys(iprow,j) / tmp
         enddo

c        --- remove from all other rows ---
         do i = 1,nrows
            if (i .ne. iprow) then
               tmp = rsys(i,ipcol)/rsys(iprow,ipcol)
               do j = 1,ncols
                  rsys(i,j) = rsys(i,j) - tmp*rsys(iprow,j)
               enddo
            endif
         enddo

c      write(6,*) 'iprow = ',iprow,' ipcol = ',ipcol
c      write(6,*) ' -------------------------------------- '
c      do i = 1,nrows
c         write(6,'(20f6.2)') rsys(i,:)
c      enddo

	enddo
 201  continue
c      write(6,*) ' ****** END PHASE I ***'

c     ----- perform phase II simplex -----
      do k = 1,10*ncols
         ipcol = 0
         iprow = 0
         tmp = 0.0d0

c        --- pivot column is most negative cost function ---
         do j = 1,2*nspecies+nactive
            if (rsys(nrows-1,j) .lt. tmp) then
               tmp = rsys(nrows-1,j)
               ipcol = j
            endif
         enddo

         if (ipcol .eq. 0) goto 202

c        --- pivot row has the smallest positive ratio ---
         tmp = 0.0d0
         do j = 1,nrows - 2
            if ( (rsys(j,ipcol) .gt. 0.0d0) .and. 
     &           (rsys(j,ncols) .ge. 0.0d0) ) then
               if ((tmp .gt. (rsys(j,ncols) / rsys(j,ipcol))) .or.
     &              (iprow .eq. 0)) then
                  tmp = rsys(j,ncols) / rsys(j,ipcol)
                  iprow = j
               endif
            endif
         enddo
         
         if (iprow .eq. 0) goto 202

c        --- normalize pivot row ---
         tmp = rsys(iprow,ipcol)
         do j = 1,ncols
            rsys(iprow,j) = rsys(iprow,j) / tmp
         enddo

c        --- remove from all other rows ---
         do i = 1,nrows-1
            if (i .ne. iprow) then
               tmp = rsys(i,ipcol)/rsys(iprow,ipcol)
               do j = 1,ncols
                  rsys(i,j) = rsys(i,j) - tmp*rsys(iprow,j)
               enddo
            endif
         enddo

c      write(6,*) 'iprow = ',iprow,' ipcol = ',ipcol
c      write(6,*) ' -------------------------------------- '
c      do i = 1,nrows
c         write(6,'(20f6.2)') rsys(i,:)
c      enddo

      enddo
 202  continue      

c     ----- save basis species -----
      do j = 1,nspecies
         if (rsys(nrows-1,j) .eq. 0.0d0) then
            k = 0
            do i = 1,nelements
               if (isactive(i) .gt. 0) then
                  k = k + 1
                  if (rsys(k,j) .gt. 0.0d0) then
                     ibasis(i) = j
                     goto 203
                  endif   !if element = 1.0 (>0.0)
               endif   !if isactive
            enddo   !do i loop
         endif  !cost function = 0.0
 203     continue
c         if (rsys(nrows-1,j) .eq. 0.0d0)
c     &                write(6,*) 'SPECIES ',j,' BASIC!'
      enddo !do j loop

      if (allocated(rsys)) deallocate(rsys)

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine simplex
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ---------------Subroutine generate_element_table()-------------------
c ---------------------------------------------------------------------
      Subroutine generate_element_table(nspecies, snames, etable, 
     &                                  nelements, ierror)
c
c *********************************************************************
c ******    generate the elemental breakdown table               ******
c ******                                                         ******
c ****** written by:    Matthew MacLean                          ******
c ******                maclean@cubrc.org                        ******
c ****** coded on:      June     13, 2009                        ******
c ****** last updated:  April    15, 2010                        ******
c *********************************************************************
c
      Implicit none

c     ----- declaration of external arguments -----
      Integer :: nspecies, nelements, ierror
      Character*10,dimension(nspecies)   :: snames
      Real*8,dimension(:,:),pointer        :: etable

c     ----- declaration of internal variables -----
      Integer :: i,j,is,ic,itype,ntmpelem
      Integer :: iclose, ipopen, iqstart, iqend, isurf,
     &           istart, iend, ival, ncurrent, ncharge, ilen
      Character :: c
      Character*10,allocatable,dimension(:) :: tmp_elem_names
      Integer,allocatable,dimension(:,:)    :: tmp_elem_matrix

c     ----- declaration of external function types -----
      Integer :: isupper, getchartype, matchelemname, isplusminus

c 
c **************** BEGIN SUBROUTINE BODY *******************
c

c     ----- initialize error to zero -----
      ierror = 0

c     ----- precompute the number of possible species that may occur -----
c           n = 1 to start to account for electrons
c           anything with 's' might be a surface site
      ntmpelem = 1
      do i = 1,nspecies
         do j = 1,10
            c = snames(i)(j:j)
            if (c .eq. 's') ntmpelem = ntmpelem + 1
            if (isupper(c) .eq. 1) ntmpelem = ntmpelem + 1
         enddo
      enddo


c     ----- declare some temporary arrays to track elements with -----
      allocate(tmp_elem_names(ntmpelem),
     &         tmp_elem_matrix(nspecies,ntmpelem),stat=ierror)
      tmp_elem_matrix = 0

c     ----- initialize the number of actual unique elements -----
      nelements = 0
      
      iclose  = 0
      ipopen  = 0
      istart  = 0
      iend    = 0
      iqstart = 0
      iqend   = 0
      isurf   = 0
      
      do is = 1,nspecies
         do ic = 1,10

c          ^^^ get the type of character we look at ^^^
           itype = getchartype(snames(is)(ic:ic))

       
           if (ipopen .gt. 0) then
c          ^^^ special read until we close parentheses ^^^    
              if (itype .eq. 7) then
                 if (ic .gt. (ipopen + 1)) then
                    if (snames(is)(ipopen+1:ipopen+1) .ne. 'b') then !don't count bulk "b" as an element
                  ncurrent = matchelemname(nelements, tmp_elem_names,
     &                        snames(is)(ipopen+1:ic-1), ic-ipopen-1)
                       tmp_elem_matrix(is,ncurrent) = 1    !default stoichiometry
			                 endif
                    ipopen = 0
                 else
                    ipopen = 0  !just close - empty parentheses
                 endif
              endif
           else
c          ^^^ case statement to parse string ^^^    
           select case (itype)
           
           case (1) !Upper Case - starts new element
              if (istart .gt. 0) then   !already an element
                 ncurrent = matchelemname(nelements, tmp_elem_names, 
     &                         snames(is)(istart:iend), iend-istart+1)
                 tmp_elem_matrix(is,ncurrent) = 1    !default stoichiometry
              endif
              
              if (snames(is)(ic:ic) .ne. 'E') then
                istart = ic
                iend   = ic
              else
                istart = 0
              endif
              
           case (2) !lower Case - continue
           
              if (istart .gt. 0) then
                 iend   = iend + 1
              else
                 if (snames(is)(ic:ic) .ne. 'e') then !free electron
                    ncharge = matchelemname(nelements, tmp_elem_names, 
     &                         'CHARGE', 6)
                    tmp_elem_matrix(is,ncharge) = 
     &                                tmp_elem_matrix(is,ncharge) - 1
                 else  !error - no other elements start with lowercase
                   ierror = 2
                 endif
              endif
           
           case (3) !Number - quantity info
           
              if (iqstart .eq. 0) then
                 iqstart = ic
                 iqend = ic
              else
                 iqend = iqend + 1
              endif
           
           case (4:5) !Plus or Minus - charge
              
              if (istart .gt. 0) then   !already an element
                 ncurrent = matchelemname(nelements, tmp_elem_names, 
     &                         snames(is)(istart:iend), iend-istart+1)
                 tmp_elem_matrix(is,ncurrent) = 1    !default stoichiometry
                 istart = 0
              endif
              
              ncharge = matchelemname(nelements, tmp_elem_names, 
     &                         'CHARGE', 6)
              tmp_elem_matrix(is,ncharge) = 
     &                                tmp_elem_matrix(is,ncharge) 
     &                              + isplusminus(snames(is)(ic:ic))
           
           case (6)
              
              if (ipopen .gt. 0) then
                 ierror = 6
              else
                 ipopen = ic
              endif

           case (7)

              if (ipopen .eq. 0) then
                 ierror = 7
              else
                 ipopen = 0
              endif
           
           case default
           
           end select
c          ^^^ finished select ^^^

c          ^^^ determine if we need to read a species ^^^
           if ((istart .gt. 0).and.
     &         ((itype .ne. 1).and.(itype .ne. 2)) ) iclose = 1
           if ((istart .gt. 0).and.(ic .eq. 10)) iclose = 1

           if (iclose .eq. 1) then
              
              ncurrent = matchelemname(nelements, tmp_elem_names, 
     &                         snames(is)(istart:iend), iend-istart+1)
              tmp_elem_matrix(is,ncurrent) = 1    !default stoichiometry
              iclose = 0
              istart = 0
           endif

c          ^^^ determine if we need to read a quantity ^^^
           if ((iqstart .gt. 0).and.(itype .ne. 3)) iclose = 1
           if ((iqstart .gt. 0).and.(ic .eq. 10)) iclose = 1

           if (iclose .eq. 1) then
              read(snames(is)(iqstart:iqend),'(i10)') ival
              tmp_elem_matrix(is,ncurrent) = ival
              iclose = 0
              iqstart = 0
           endif

           endif

         enddo
      enddo    

      write(6,*) '     Element List '
      write(6,'(a,$)') '       '
      do i = 1,nelements
         ilen = len_trim(tmp_elem_names(i))
         write(6,'(a,a,$)') '  ',tmp_elem_names(i)(1:ilen)
      enddo
      write(6,*)
      
c     ----- copy over to permanent arrays -----
      allocate(etable(nspecies,nelements),stat=ierror)
      etable = 0
      do i = 1,nspecies
         do j = 1,nelements
            etable(i,j) = dble(tmp_elem_matrix(i,j))
         enddo
      enddo

      write(6,*) '     Stoichiometry'
      do i = 1,nspecies
         write(6,200) snames(i),etable(i,1:nelements)
      enddo
 200  format('         ',a,':',20f5.0)

      if (allocated(tmp_elem_names))  deallocate (tmp_elem_names)
      if (allocated(tmp_elem_matrix)) deallocate (tmp_elem_matrix)

c 
c **************** END SUBROUTINE BODY *********************
c
      End Subroutine generate_element_table
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ----------------Function getchartype()---------------------
c ---------------------------------------------------------------------
c     RETURN TYPE TABLE
c     0 = Unrecognized
c     1 = Upper Case Character
c     2 = Lower Case Character
c     3 = Number
c     4 = Plus Sign "+"
c     5 = Minus Sign "-"
c     6 = Left Parenthesis "("
c     7 = Right Parenthesis ")"
c     9 = Space " "
c

      Function getchartype(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, getchartype

      iascii = ichar(c)
      
      if ((iascii .ge. 65).and.(iascii .le. 90)) then
         getchartype = 1
      elseif ((iascii .ge. 97).and.(iascii .le. 122)) then
         getchartype = 2
      elseif ((iascii .ge. 48).and.(iascii .le. 57)) then
         getchartype = 3
      elseif (iascii .eq. 43) then
         getchartype = 4
      elseif (iascii .eq. 45) then
         getchartype = 5
      elseif (iascii .eq. 40) then
         getchartype = 6
      elseif (iascii .eq. 41) then
         getchartype = 7
      elseif (iascii .eq. 32) then
         getchartype = 8
      else
         getchartype = 0
      endif
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function getchartype
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ----------------Function matchelemname()---------------------
c ---------------------------------------------------------------------
      Function matchelemname(nelements, elem_names, cname, nc)

      Implicit none
      
      Character*(*) :: cname
      Integer       :: nelements, nc
      Character*10,dimension(nelements+1) :: elem_names
      Integer     :: i, matchelemname

c     ----- return an index indicating a new element by default -----
      matchelemname = nelements + 1

      do i = 1,nelements
        if (elem_names(i)(1:nc) .eq. cname(1:nc)) matchelemname = i
      enddo

      if (matchelemname .gt. nelements) then
         nelements = nelements + 1
         elem_names(nelements) = cname(1:nc)
      endif
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function matchelemname
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ----------------Function isupper()---------------------
c ---------------------------------------------------------------------
      Function isupper(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, isupper


      iascii = ichar(c)
      
      if ((iascii .ge. 65).and.(iascii .le. 90)) then
         isupper = 1
      else
         isupper = 0
      endif
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function isupper
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ----------------Function islower()---------------------
c ---------------------------------------------------------------------
      Function islower(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, islower


      iascii = ichar(c)
      
      if ((iascii .ge. 97).and.(iascii .le. 122)) then
         islower = 1
      else
         islower = 0
      endif
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function islower
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ----------------Function isnumber()---------------------
c ---------------------------------------------------------------------
      Function isnumber(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, isnumber


      iascii = ichar(c)
      
      if ((iascii .ge. 48).and.(iascii .le. 57)) then
         isnumber = iascii - 48
      else
         isnumber = -1
      endif
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function isnumber
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-




c --------------------------------------------------------------------- 
c ----------------Function isparenthesis()---------------------
c ---------------------------------------------------------------------
      Function isparenthesis(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, isparenthesis


      iascii = ichar(c)
      
      select case (iascii)
      case(40)
         isparenthesis = 1
      case(41)
         isparenthesis = -1
      case default
         isparenthesis = 0
      end select
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function isparenthesis
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



c --------------------------------------------------------------------- 
c ----------------Function isparenthesis()---------------------
c ---------------------------------------------------------------------
      Function isplusminus(c)

      Implicit none
      
      Character*1 :: c
      Integer     :: iascii, isplusminus


      iascii = ichar(c)
      
      select case (iascii)
      case(43)
         isplusminus = 1
      case(45)
         isplusminus = -1
      case default
         isplusminus = 0
      end select
c 
c **************** END SUBROUTINE BODY *********************
c
      End Function isplusminus
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
