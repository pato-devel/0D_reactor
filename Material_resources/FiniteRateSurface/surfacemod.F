      module finiteRateSurface

c     ============== Array sizing parameters (hardwired) ===============
      Integer :: ISMODEL_MAX_REACT, ISMODEL_MAX_PROD
      parameter(ISMODEL_MAX_REACT=3,ISMODEL_MAX_PROD=3)
      Integer :: NLEWIS_FITS, NLEWIS_COEFFS
      parameter(NLEWIS_FITS=6, NLEWIS_COEFFS=10)

c     ================ Version parameters (hardwired) ==================
      Integer :: MAJOR_FRM_VERSION, MINOR_FRM_VERSION
      parameter(MAJOR_FRM_VERSION=3,MINOR_FRM_VERSION=9)

c     ======= Defined Type to store constants for a surface model ======
      type finiteRateSurfaceModelType

c       ---- title
        character*200 :: title

c       ---- initialization flag(s)
        integer :: initsurf

c       ---- phase and species information
        integer ngp,nsp,nbp,nptot,kstot,ngps,nphases,nspecies,nspasmax,
     &          igstart,igend,isstart,isend,ibstart,ibend
        integer,pointer,dimension(:)   :: nsps,nbps,ithermo,nspas,
     &                                    kphase,kthermo,kadsr,ksites
        integer,pointer,dimension(:,:) :: nspass

c       ---- Thermodynamic Lewis data ----
        real*8,pointer,dimension(:,:,:) :: rlewis_coeffs,rlewis_trange
        integer,pointer,dimension(:)    :: rlewis_nrange

c       ---- Reaction information
        integer :: isrt
        integer,pointer,dimension(:,:) :: RR,RP,RVR,RVP
        integer,pointer,dimension(:,:) :: vr,vp,v
        integer,pointer,dimension(:)   :: isrfon,isrbon,istype,
     &                                    iadtype,iadform

c       ---- Phase and species information
        real*8 :: sdentot,sfrcI,vfbI,bdentot,bportot
        real*8,pointer,dimension(:)   :: sden,sfrc,den,por,vfb
        real*8,pointer,dimension(:,:) :: sumv,bmf,sdenas
        real*8,pointer,dimension(:)   :: wmass,Ediss,Ed,Eact

c       ---- Equation information
        real*8,pointer,dimension(:)   :: Cf,eta,Ea,sitep
        real*8,pointer,dimension(:)   :: Cfad,etaad,vad,Edes
        real*8,pointer,dimension(:)   :: sdensity,sfraction

c       ---- Pyrolysis blowing information
        integer :: nblwflag,nebc,nblw,isblwon
        real*8  :: cyield,ha
        real*8,pointer,dimension(:)   :: rmdotg,pyroadd
        real*8,pointer,dimension(:,:) :: blwmf

c       ---- Species and phase names
        character*10,pointer,dimension(:)   ::  species,phases

c       ---- debugging flag ---- 
        integer :: ifrm_debug 

c       ---- error string ----
        character*120 :: cfrm_error_string

      end type


c     ======== Defined Type to store concentrations for 1 cell =========
      type finiteRateSurfaceCell

c        ---- pointer to the model record 
         type(finiteRateSurfaceModelType),pointer :: frs_record

c        ---- parameters and governing fields
         real*8  :: gas_tmpr, wall_tmpr, wall_pres
         Integer :: num_species, num_gas, num_surf, num_nongas,
     &              num_nonbulk, num_eqns

c        ---- matrices for system
         real*8,pointer,dimension(:)   :: concentration, delta_conc
         real*8,pointer,dimension(:)   :: srhs
         real*8,pointer,dimension(:,:) :: sjacobian

      end type


      end modulefiniteRateSurface
