      program srmodel
c
c *********************************************************************
c ******   a standalone code to compute finite-rate gas/surface  ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      May       9, 2009                        ******
c ****** last updated:  September 6, 2011                        ******
c *********************************************************************
c      
      Use finiteRateSurface
      Implicit none

      real*8 :: PI, R, AV, PREF
      parameter(PI=3.1415926536,R=8.31451,AV=6.0221367e23,PREF=1.0e5)
	integer :: IRELEASE
	parameter(IRELEASE=12)
      
      integer :: ierr, i, k, ngas, nsteps, istep, itorder, idiagnostic,
     &           iscan, iout, iconv, igastype, ieqtype, nsphase,
     &           iscanflag, j, n, iscanextra
      real*8  :: tw, tg, p, dt, time, dti, tmp, resL1, resL2, alpha,
     &           beta, gamma, smoltot, v, volume, new_volume, sblow,
     &           pblow, final_pres, resid, sscanp, sscantw, sscantg
      character*10  :: reactnum
	character*40  :: pname, gname
	character*128 :: fname, blwname
	character*200 :: bname
      character*6,allocatable,dimension(:) :: snames
      real*8,allocatable,dimension(:) :: smolefr, eta, coftmpr, dxold,
     &                                   sscanmole
      real*8,allocatable,dimension(:,:,:) :: cof
	
	Integer :: finite_rate_surf_sanity_check
      Real*8  :: finite_rate_surf_solve_qss
      
      type(finiteRateSurfaceModelType) :: frmr
      type(finiteRateSurfaceCell)      :: cell, qss



c-------------------------------------------------
c-----Write Title and Version Info to Screen
c-------------------------------------------------

      write(6,100) 'SRModel'
      write(6,101) MAJOR_FRM_VERSION, MINOR_FRM_VERSION, IRELEASE
      write(6,102) 'written by: ','Joe Marschall  ',
     &             'jochen.marschall@sri.com'
      write(6,102) '','Matthew MacLean','maclean@cubrc.org'
      write(6,*) 



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
      read (13,*) tw, p, igastype, iscan
      tg = tw  !we no longer allow separate wall and gas temperatures

      read (13,*)
      read (13,*)
      read (13,*) dt, nsteps, itorder, idiagnostic, ieqtype

      read (13,*)
      read (13,*)
      read (13,*) ngas

      allocate(snames(ngas),smolefr(ngas),sscanmole(ngas))

      read (13,*)
      read (13,*)
      do i = 1,ngas
         read (13,*) snames(i), smolefr(i)
      enddo
      read (13,*)

c     ----- limit flags -----
      if ((itorder .ne. 2).and.(itorder .ne. 0)) itorder = 1

      if ((idiagnostic .gt. 3).or.(idiagnostic .lt. 1)) idiagnostic=0

	if (igastype .lt. 0) igastype = 0
	if (igastype .gt. 2) igastype = 2

      if ((ieqtype .ge. 3).or.(ieqtype .lt. 0)) ieqtype = igastype

c     ----- set scanning mode (turns off many flags) -----
      if (iscan .gt. 0) then
         iscan = 1        !start with scan #1
         iscanextra = 0   !turn off extra reaction output
         if (idiagnostic .gt. 0) iscanextra = 1
         idiagnostic = 0  !turn off all diagnostic output
         iout = -1        !turn off all file output
         iconv = -1       !turn off convergence output
         iscanflag = 1
      endif



c-------------------------------------------------
c-----Read in input files and collect data
c-------------------------------------------------

c     ----- read in data from surface file -----
      call read_finite_rate_surf_file(fname, frmr, ierr)
      if (ierr .gt. 0) then
         write(6,'(a,a,i3,a)')
     &           'ERROR: Call to read_finite_rate_surf_file',
     &           ' failed with error code ',ierr,'. Execution halted.'
         stop
      endif
      if (frmr%ngps .ne. ngas) then
         write(6,*) 'ERROR: Number of species in input file ',
     &              'must match chemistry database!'
         stop
      endif

c     ----- read pyrolysis blowing data if necessary -----
      if (frmr%nblwflag .gt. 0) then
c        look for a . extension on the model file
         k = scan(fname,'.',.TRUE.)
         if (k .eq. 0) then
            k = LEN_TRIM(blwname)
         else
            k = k - 1 !remove the .
         endif
c        look for the path of the model file
         i = scan(fname,'/',.TRUE.)
         i = i + 1   !go past the / or move to first character
c        try stripping the file to look in the "current" directory first
         blwname = fname(i:k)//'.blw'
         call read_finite_rate_blow_file(blwname, frmr, ierr)
c        try reading the file with the full path of the input
         if (ierr .gt. 0) then
            blwname = fname(1:k)//'.blw'
            call read_finite_rate_blow_file(blwname, frmr, ierr)
         endif
c        didn't find it -- I give up!
         if (ierr .gt. 0) then
            write(6,*) 'WARNING: Blowing data file could not be read.',
     &                 ' Pyrolysis will be de-activated.'
         endif
      endif

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
      call read_lewis_bulk_file(bname, frmr, cof, coftmpr, 5, ierr)
	if (ierr .gt. 0) then
         write(6,'(a,a,i3,a)')
     &           'ERROR: Call to read_lewis_bulk_file',
     &           ' failed with error code ',ierr,'. Execution halted.'
         stop
      endif

c     ----- write the model back out as a sanity check -----
	if (idiagnostic .ge. 3) then
         open(10,file='model.out')
         call write_finite_rate_surf(frmr, 10, ierr)
         close(10)
      endif

c     ----- perform sanity checks on the data -----
      ierr = finite_rate_surf_sanity_check(frmr)
      if ((ierr .gt. 0).and.(ierr .lt. 100)) 
     &   write(6,103) ierr,trim(frmr%cfrm_error_string)
      if (ierr .ge. 100) then
         write(6,104) ierr,trim(frmr%cfrm_error_string)
         stop
      endif

c     ----- open an output file if needed -----
      if (iout  .lt. 0) iout = 0
      if (iconv .lt. 0) iout = 0
      if (iout  .ne. 0) open(14, file=pname)


c-------------------------------------------------
c-----Initialize single surface cell
c-------------------------------------------------

c     ----- initialize and allocate arrays -----
      call finite_rate_surf_initialize(cell, frmr)

c     ----- set values by hand -----
      cell%gas_tmpr  = tg
      cell%wall_tmpr = tw
      cell%wall_pres = p

c     ----- initialize unit volume -----
      volume = 1.0d0
      new_volume = volume

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

c     ----- initialize to QSS if needed -----
      if (frmr%initsurf .eq. 1) resid = finite_rate_surf_solve_qss(cell)

c     ----- write initial cell state -----
      write(6,*)
      write(6,110)
c      write(6,111) 'Gas Tmpr   ',cell%gas_tmpr,'K'
      write(6,111) 'Wall Tmpr  ',cell%wall_tmpr,'K'
      write(6,111) 'Pressure   ',cell%wall_pres,'Pa'
      do k = 1,cell%num_gas
         write(6,111) 'x'//frmr%species(k),
     &                cell%concentration(k),'mol/m^3'
      enddo
      do k = cell%num_gas+1,cell%num_nonbulk
         write(6,111) 'x'//frmr%species(k),
     &                cell%concentration(k),'mol/m^2'
      enddo
      do k = cell%num_nonbulk+1,cell%num_species
         write(6,111) 'x'//frmr%species(k),
     &                cell%concentration(k),'[]'
      enddo

c     ----- finally save a copy of initial values for equilibrium calc -----
      allocate(eta(cell%num_species))
      do k=1,cell%num_species
         eta(k) = cell%concentration(k)
      enddo



c-------------------------------------------------
c-----Set time integration flags
c-------------------------------------------------

      if (itorder .eq. 2) then
         alpha = 1.5d0
         beta  = 0.5d0
         gamma = 1.0d0
      elseif (itorder .eq. 1) then
         alpha = 1.0d0
         beta  = 0.0d0
         gamma = 1.0d0
      else
         alpha = 1.0d0
         beta  = 0.0d0
         gamma = 0.0d0
      endif
      
      dti  = alpha / dt
      beta = beta / dt

      allocate(dxold(cell%num_eqns))
      dxold = 0.0d0



c-------------------------------------------------
c-----Write header info and track scan mode
c-------------------------------------------------

      write(6,*)
      if (iscan .le. 0) then
         write(6,113)
         write(6,201)
      endif

      if (iout .gt. 0) then
         write(14,301)
         do k = 1,cell%num_species
            i = len_trim(frmr%species(k))
            write(14,302) 'x'//frmr%species(k)(1:i)
         enddo
         write(14,*)
         write(14,310) 0, 0.0
         do k = 1,cell%num_species
            write(14,311) cell%concentration(k)
         enddo
         write(14,*)
      endif
      
 901  continue    !return here on scan mode only & reset cell
      if (iscan .gt. 0) then

c        ----- reset time integration terms -----
         dxold      = 0.0d0
         volume     = 1.0d0
         new_volume = volume

c        ----- read new line of settings -----
         if (iscanflag .eq. 1) then
           read(13,*)
           read(13,*)
           read(13,*)
           iscan = iscan + 1
           iscanflag = 0
           select case (iscan)
           case(2)
              write(6,501)
              write(6,508) 'Pi','Tw','Volume','Pf','resL2'
           case(3)
              write(6,502)
              do k = 1,ngas
                i = len_trim(frmr%species(k))
                write(6,509) 'mfi_'//frmr%species(k)(1:i)
              enddo
              write(6,509) 'Volume'
              write(6,509) 'Pf'
		  write(6,509) 'resL2'
           end select
           if (iscan .le. 3) then
             do k = 1,cell%num_species
               i = len_trim(frmr%species(k))
               write(6,509) 'x'//frmr%species(k)(1:i)
             enddo
             do k = 1,cell%num_gas
               i = len_trim(frmr%species(k))
               write(6,509) 'g'//frmr%species(k)(1:i)
             enddo
             if (iscanextra .gt. 0) then
             do k = 1,cell%num_species
               i = len_trim(frmr%species(k))
               do j = 1,frmr%isrt
                 write(reactnum,'(i10)') j
                 reactnum = adjustl(reactnum)
                 n = len_trim(reactnum)
                 write(6,509)
     &                 'w'//frmr%species(k)(1:i)//'_'//reactnum(1:n)
               enddo
             enddo
             endif
           endif
           write(6,*)
           
         endif

c        ----- read a new row of data -----
         select case (iscan)
         case (2)   !p, tg, tw
            read(13,*) sscanp, sscantw
            sscantg = sscantw  !no longer allow separate wall and gas temperatures
            if ((sscantg .lt. 0.0d0) .or. (sscantw .lt. 0.0d0)
     &                            .or. (sscanp .lt. 0.0d0)) then
               iscanflag = 1
               go to 901
            endif
            cell%gas_tmpr  = sscantg
            cell%wall_tmpr = sscantw
            cell%wall_pres = sscanp
            sscanmole      = smolefr  !copy default mole fractions
         case (3)   !mole fractions
            read(13,*) sscanmole(1:ngas)
            smoltot = 0.0d0
            do i = 1,ngas
              if (sscanmole(i) .lt. 0.0d0) then
                 iscanflag = 1
                 go to 901
              endif
              if (sscanmole(i) .lt. 1.0d-15) sscanmole(i) = 1.0d-15
              smoltot = smoltot + sscanmole(i)
            enddo
            do i = 1,ngas
              sscanmole(i) = sscanmole(i)/smoltot
            enddo
            sscantg = tg
            sscantw = tw
            sscanp  = p
            cell%gas_tmpr  = tg  !copy defaults for tmpr and pressure
            cell%wall_tmpr = tw
            cell%wall_pres = p
         case default
            go to 902
         end select

c        ----- reset gas cell concentrations -----
         do k=1,ngas
           cell%concentration(k) = cell%wall_pres *
     &                                sscanmole(k)/(R*cell%gas_tmpr)
         enddo
         do k=ngas+1,cell%num_species  !surface concentrations remain the same
            cell%concentration(k) = eta(k)
         enddo
c        ----- initialize to QSS if needed -----
         if (frmr%initsurf .eq. 1)
     &       resid = finite_rate_surf_solve_qss(cell)
      endif



c-------------------------------------------------
c-----Integrate in time
c-------------------------------------------------

      time = 0.0d0

      do istep = 1,nsteps

c        ----- compute implicit system -----
         call finite_rate_compute_system(cell)

c        ----- freeze gas phase at user request -----
         if (igastype .eq. 1) then
            do k = 1,cell%num_gas
               do i = 1,cell%num_eqns
                  cell%sjacobian(k,i) = 0.0d0
               enddo
               cell%sjacobian(k,k) = 1.0d0
               cell%srhs(k) = 0.0d0
            enddo
         endif

c        ----- freeze bulk phase (always) -----
         do k = cell%num_nonbulk+1,cell%num_species
            do i = 1,cell%num_eqns
               cell%sjacobian(k,i) = 0.0d0
            enddo
            cell%sjacobian(k,k) = 1.0d0
            cell%srhs(k) = 0.0d0
         enddo

c        ----- add time step term and negate -----
         do k = 1,cell%num_eqns
            do i = 1,cell%num_eqns
               if (i .eq. k) then
                  cell%sjacobian(i,k) = dti - gamma*cell%sjacobian(i,k)
               else
                  cell%sjacobian(i,k) = - gamma*cell%sjacobian(i,k)
               endif
            enddo
         enddo

c        ----- add 2nd order time explicit term -----
         do k = 1,cell%num_eqns
            cell%srhs(k) = cell%srhs(k) + beta * dxold(k)
         enddo

c        ----- solve system -----
         call lusolve(cell%sjacobian, cell%srhs, cell%num_eqns)

c        ----- compute new volume for constant pressure -----
         if (igastype .eq. 2) then
           new_volume = 0.0d0
           do k = 1,cell%num_gas
              tmp = cell%srhs(k) + volume * cell%concentration(k)
              if (tmp .lt. 1.0d-20) tmp = 1.0d-20
              new_volume = new_volume + tmp
           enddo
           new_volume = new_volume * (R*cell%gas_tmpr/cell%wall_pres)
         else
           new_volume = volume
         endif

c        ----- update variables (accounting for volume change) -----
         do k = 1,cell%num_gas
            tmp = volume * cell%concentration(k)
            cell%concentration(k) = (tmp + cell%srhs(k))/new_volume
            if (cell%concentration(k) .lt. 1.0d-20)
     &             cell%concentration(k) = 1.0d-20
         enddo
         do k = cell%num_gas+1,cell%num_species
            cell%concentration(k) = cell%concentration(k) + cell%srhs(k)
            if (cell%concentration(k) .lt. 1.0d-20)
     &             cell%concentration(k) = 1.0d-20
         enddo
         cell%gas_tmpr  = cell%gas_tmpr  + cell%srhs(cell%num_species+1)
         cell%wall_tmpr = cell%wall_tmpr + cell%srhs(cell%num_species+2)
         volume         = new_volume

c        ----- increment time index -----
         time = time + dt

c        ----- save solution for next iteration (2nd order time) -----
         dxold = cell%srhs

c        ----- compute residual -----
         resL1 = 0.0d0
         resL2 = 0.0d0
         do k = 1,cell%num_species
            tmp = abs(cell%srhs(k))
            if (tmp .gt. resL1) resL1 = tmp
            resL2 = resL2 + tmp * tmp
         enddo
         resL2 = sqrt(resL2)

c        ----- write output to screen & file -----
         if (iconv .gt. 0) then
            if ((mod(istep,iconv) .eq. 0).or.(istep.eq.nsteps))
     &                    write(6,202) istep,time,resL1,resL2
         else
            if ((istep.eq.nsteps).and.(iscan .le. 0)) 
     &                    write(6,202) istep,time,resL1,resL2
         endif
         
         if(iout .gt. 0) then
            if (mod(istep,iout) .eq. 0) then
               write(14,310) istep,time
               do k = 1,cell%num_species
                  write(14,311) cell%concentration(k)
               enddo
               write(14,*)
            endif
         endif

      enddo

c     ----- compute final gas pressure -----
      final_pres = 0.0d0
      do k=1,ngas
        final_pres = final_pres + cell%concentration(k)*R*cell%gas_tmpr
      enddo

c     ----- summarize final conditions to screen -----
      if (iscan .le. 0) then
         write(6,*)
         write(6,112)
c         write(6,111) 'Gas Tmpr   ',cell%gas_tmpr,'K'
         write(6,111) 'Wall Tmpr  ',cell%wall_tmpr,'K'
         write(6,111) 'Pressure   ',final_pres,'Pa'
         write(6,111) 'Rel Volume ',volume,''
         do k = 1,cell%num_gas
            write(6,111) 'x'//frmr%species(k),
     &                   cell%concentration(k),'mol/m^3'
         enddo
         do k = cell%num_gas+1,cell%num_nonbulk
            write(6,111) 'x'//frmr%species(k),
     &                   cell%concentration(k),'mol/m^2'
         enddo
         do k = cell%num_nonbulk+1,cell%num_species
            write(6,111) 'x'//frmr%species(k),
     &                   cell%concentration(k),'[]'
         enddo
      else
         select case (iscan)
         case(2)
           write(6,510) sscanp,sscantw,volume,final_pres,resL2
         case(3)
           do k = 1,ngas
             write(6,519) sscanmole(k)
           enddo
           write(6,519) volume
           write(6,519) final_pres
	     write(6,519) resL2
         end select
         do k = 1,cell%num_species
            write(6,519) cell%concentration(k)
         enddo
         call finite_rate_compute_system(cell)
         do k = 1,cell%num_gas
            v = sqrt(R*cell%gas_tmpr/(2.0d0*PI*frmr%wmass(k)))
            write(6,519) -cell%srhs(k)/cell%concentration(k)/v
         enddo
         if (iscanextra .gt. 0) call srmodel_scan_diagnostics(cell)
         write(6,*)
         go to 901
      endif

 902  continue   !iscan comes here when done

      if (iout .ne. 0) close(14)



c-----------------------------------------------------
c-----Post-processed quantities
c-----------------------------------------------------

      if (iscan .le. 0) then

c        ----- calculate loss rates (recombination probability) -----
         write(6,*)
         write(6,114)
         call finite_rate_compute_system(cell)
         do k = 1,cell%num_gas
            v = sqrt(R*cell%gas_tmpr/(2.0d0*PI*frmr%wmass(k)))
            write(6,111) 'g'//frmr%species(k),
     &                   -cell%srhs(k)/cell%concentration(k)/v
         enddo

c        ----- calculate blowing and recession rates -----      
         write(6,*)
         write(6,115)
         sblow = 0.0d0
         pblow = 0.0d0
         do k = cell%num_nonbulk+1,cell%num_species
            sblow = sblow - frmr%wmass(k) * cell%srhs(k)
         enddo
         if (frmr%nblwflag .eq. 1) then
            if (frmr%nebc .gt. 0) then
               pblow = sblow * (1.0d0/frmr%cyield - 1.0d0)
            else
               do k = 1,cell%num_gas
                  pblow = pblow + frmr%pyroadd(k)
               enddo
            endif
         endif
         write(6,111) 'blowing rate   ',sblow+pblow,' kg/m^2-s'
         write(6,111) 'recession rate ',
     &                    sblow/max(frmr%bdentot,1.0d-30),' m/s'

      endif



c-------------------------------------------------
c-----Perform QSS (Quasi-Steady State) Analysis
c-------------------------------------------------

      if (idiagnostic .ge. 1) then

c        ----- initialize and allocate arrays -----
         call finite_rate_surf_initialize(qss, frmr)

c        ----- write header -----
         write(6,*) 
         write(6,130)

c        ----- set pressure, temperatures, and gas species -----
         qss%gas_tmpr  = tg
         qss%wall_tmpr = tw
         qss%wall_pres = p
         do k = 1,cell%num_gas
            qss%concentration(k) = cell%concentration(k)
         enddo

c        ----- call QSS solver -----
         resid = finite_rate_surf_solve_qss(qss)

c        ----- write to screen -----
         write(6,131) sqrt(resid)
         write(6,132)
         do k = cell%num_gas+1,cell%num_nonbulk
            write(6,133) 'x'//frmr%species(k),qss%concentration(k)
         enddo

      endif



c-------------------------------------------------
c-----Perform Equilibrium Computation
c-------------------------------------------------

      if (idiagnostic .ge. 1) then

c        ----- write header and call solver -----
         write(6,*) 
         write(6,120)
         call srmodel_compute_equilibrium(cell, eta, ieqtype)
         write(6,121)

c        ----- write concentrations -----      
         do k = 1,cell%num_gas
            write(6,122) 'x'//frmr%species(k),eta(k),'mol/m^3'
         enddo
	   do k = cell%num_gas+1,cell%num_nonbulk
            write(6,122) 'x'//frmr%species(k),eta(k),'mol/m^2'
         enddo
	   do k = cell%num_nonbulk+1,cell%num_species
            write(6,122) 'x'//frmr%species(k),eta(k),'[]'
         enddo

      endif



c-------------------------------------------------
c-----Perform Extra Convergence Analysis
c-------------------------------------------------

      if (idiagnostic .ge. 2) call srmodel_compute_diagnostics(cell)



c-------------------------------------------------
c-----Clean-up Memory & Close Input
c-------------------------------------------------

      deallocate(dxold, smolefr, coftmpr, cof, snames, eta, sscanmole)

      close(13)



c     ----- write string -----
 100  format(a)
c     ----- write version -----
 101  format('   version:     ',i1,'.0',i1,'.',i2)
c     ----- author field -----
 102  format('   ',a12,' ',a15,' (',a,')')
c     -----  Warning Codes -----
 103  format('WARNING CODE ',i3,': ',a)
c     ----- Error Codes -----
 104  format('ERROR CODE ',i3,': ',a)
c     ----- initial & final conditions -----
 110  format('---Initial Conditions-----------------------------')
 111  format('      ',a,' = ',es13.4,' ',a)
 112  format('---Final Conditions-------------------------------')
 113  format('--------------------------------------------------')
 114  format('---Loss Rates-------------------------------------')
 115  format('---Blowing & Recession Rate-----------------------')
c     ----- equilibrium solver -----
 120  format('---Gibbs Free-energy Minimization-----------------')
 121  format('      Concentrations')
 122  format('         ',a,' = ',es13.4,' ',a)
c     ----- QSS solver -----
 130  format('---QSS Surface Coverage Distribution--------------')
 131  format('      Residual = ',es13.4)
 132  format('      Concentrations')
 133  format('         ',a,' = ',es13.4,' mol/m^2')
c     ----- residual header output -----
 201  format('iteration         time       L1-norm       L2-norm')
c     ----- residual line output -----
 202  format(i8,x,es13.4,x,es13.4,x,es13.4)
c     ----- output file header -----
 301  format('iteration          time',$)
 302  format(x,a15,$)
c     ----- output file line -----
 310  format(i9,x,es13.4,$)
 311  format(x,es15.6,$)
c     ----- analysis output -----
 401  format('                Finite Rate         Equilibrium')
 402  format('EqC(',i3,') =   ',es14.6,$)
 403  format('       ',es14.6)
c     ----- scan mode header -----
 501  format('---Pressure/Temperature Scan----------------------')
 502  format('---Mole Fraction Scan----------------------')
 508  format(5a15,$)
 509  format(x,a14,$)
c     ----- scan mode output -----
 510  format(x,es14.6,x,es14.6,x,es14.6,x,es14.6,x,es14.6,$)
 519  format(x,es14.6,$)
      end program



c --------------------------------------------------------------------- 
c --------------Function finite_rate_surf_solve_qss()------------------
c ---------------------------------------------------------------------
      Real*8 Function finite_rate_surf_solve_qss(qss)
c
c *********************************************************************
c ******    compute QSS surface coverage holding the gas         ******
c ******    concentrations constant.                             ******
c ******    return: residual of the solver                       ******
c ******                                                         ******
c ****** written by:    Matthew MacLean (maclean@cubrc.org)      ******
c ******                Joe Marschall (jochen.marschall@sri.com) ******
c ****** coded on:      March     4, 2011                        ******
c ****** last updated:  March     4, 2011                        ******
c *********************************************************************
c

c *** provided by finiteRateSurface
c     finiteRateSurfaceCell      = type declaration for cell
c     finiteRateSurfaceModelType = type declaration for model
c
      Use finiteRateSurface
c *** no implicit types allowed
      Implicit none
c 
c *** declarations of argument types passed to the subroutine
c     qss    = initialized cell with all concentrations already set
c     smodel = copy of pointer to the surface model
      type(finiteRateSurfaceCell)   :: qss
      type(finiteRateSurfaceModelType),pointer :: smodel

c 
c **** declaration of internally used variables
c     qsslhs  = left-hand side Jacobian matrix
c     qssrhs  = right-hand side explicit matrix
c     k       = counter
c     istep   = step counter
c     nsphase = total number of active site phases
c     resid   = residual
c     rlimit  = limiter (inactive)
      Integer :: k,istep,nsphase
      Real*8  :: resid, rlimit
      Real*8,allocatable,dimension(:) :: qssrhs
      Real*8,allocatable,dimension(:,:) :: qsslhs
c 
c **** declaration of functions used
c      NONE

c 
c *************** BEGIN FUNCTION BODY ********************
c 

c     ----- initialize residual to garbage
      finite_rate_surf_solve_qss = -1.0d0

c---- copy model pointer (easy access)
      smodel => qss%frs_record

c     ----- determine number of phases -----
      nsphase = 0
      do k = 1,smodel%nsp
         nsphase = nsphase + smodel%nspas(k)
      enddo

c     ----- allocate arrays -----
      allocate(qsslhs(qss%num_surf+nsphase,qss%num_surf+nsphase),
     &         qssrhs(qss%num_surf+nsphase))

c     ----- perform Newton iteration for surface distribution -----
      do istep = 1,25
            
c        ----- get matrices -----
         call finite_rate_surf_qss_matrix(qss,qsslhs,qssrhs,
     &                                    qss%num_surf+nsphase)

c        ----- solve system -----
         call lusolve(qsslhs, qssrhs, qss%num_surf+nsphase)

c        ----- update variables -----
         resid = 0.0d0
         rlimit = 1.0d0
         do k = 1,qss%num_surf
            qss%concentration(k+qss%num_gas) = 
     &      qss%concentration(k+qss%num_gas) + rlimit*qssrhs(k)
            if (qss%concentration(k+qss%num_gas) .lt. 1.0d-20) 
     &          qss%concentration(k+qss%num_gas) = 1.0d-20
            resid = resid + qssrhs(k)*qssrhs(k)
         enddo
            
         if (resid .lt. 1.0d-30) go to 153

      enddo
 153  continue

c     ----- store residual for return -----
      finite_rate_surf_solve_qss = resid

c     ----- reclaim memory -----
      deallocate(qsslhs, qssrhs)

c 
c **************** END FUNCTION BODY *********************
c
      End Function finite_rate_surf_solve_qss
c -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-
