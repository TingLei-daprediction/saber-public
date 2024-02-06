!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_parameter
!***********************************************************************
!                                                                      !
!       Set resolution, grid and decomposition                         !
!                       - offset version -                             !
!                                                                      !
!  Note:  ixm(1)=nxm, jym(1)=mym                                       !
!                                                                      !
!    If   mod(nxm,2)=0  then    mod(im0,4)=0                           !
!    If   mod(mym,2)>0  then    mod(im0,8)=0                           !
!                                                                      !
!  Modules: kinds, jp_pietc                                            !
!                                                     M. Rancic (2022) !
!***********************************************************************
use mpi
use kinds, only: i_kind,r_kind
use jp_pietc, only: u1
!use berror, only: mg_ampl0,im_filt,jm_filt
!TEST
!use mpimod, only: nxpe,nype
!TEST

implicit none
#if 0
xxx
type mg_parameter_type
#endif
!-----------------------------------------------------------------------
!*** 
!*** Namelist parameters
!***
real(r_kind):: mg_ampl01,mg_ampl02,mg_ampl03
real(r_kind):: mg_weig1,mg_weig2,mg_weig3,mg_weig4
integer(i_kind):: mgbf_proc
logical:: mgbf_line
integer(i_kind):: nxPE,nyPE,im_filt,jm_filt
logical:: lquart

!*** 
!*** Number of generations
!***
integer(i_kind):: gm            

!*** 
!*** Horizontal resolution 
!***

!
! Original number of data on GSI analysis grid
!
integer(i_kind):: nA_max0
integer(i_kind):: mA_max0

!
! Global number of data on Analysis grid
!
integer(i_kind):: nm0        
integer(i_kind):: mm0       

!
! Number of PEs on Analysis grid
!
integer(i_kind):: nxm           
integer(i_kind):: mym           

!
! Number of data on local Analysis grid
!
integer(i_kind):: nm         
integer(i_kind):: mm        

!
! Number of data on global Filter grid
!
integer(i_kind):: im00
integer(i_kind):: jm00

!
! Number of data on local  Filter grid
!
integer(i_kind):: im
integer(i_kind):: jm    

!
! Initial index on local  Filter grid
!
integer(i_kind):: i0
integer(i_kind):: j0    
!
! Initial index on local analysis grid
!
integer(i_kind):: n0
integer(i_kind):: m0    

!
! Halo on local Filter grid 
!
integer(i_kind):: ib
integer(i_kind):: jb         

!
! Halo on local Analysis grid 
!
integer(i_kind):: nb
integer(i_kind):: mb     


integer(i_kind):: hx,hy,hz
integer(i_kind):: p
integer(i_kind):: nh,nfil
real(r_kind):: pasp01,pasp02,pasp03
real(r_kind):: pee2,rmom2_1,rmom2_2,rmom2_3,rmom2_4


integer, allocatable, dimension(:):: maxpe_fgen
integer, allocatable, dimension(:):: ixm,jym,nxy
integer, allocatable, dimension(:):: im0,jm0
integer, allocatable, dimension(:):: Fimax,Fjmax
integer, allocatable, dimension(:):: FimaxL,FjmaxL

integer(i_kind):: npes_filt

integer(i_kind):: maxpe_filt

integer(i_kind):: imL,jmL
integer(i_kind):: lm            ! number of vertical layers
integer(i_kind):: lm05          ! half of vertical levels
integer(i_kind):: km3           ! number of 3d variables
integer(i_kind):: km2           ! number of 2d variables
integer(i_kind):: km            ! number of vertically stacked all variables (km=km2+lm*km3)
integer(i_kind):: lm_all        ! vertically stacked all variables (lm_all=km - for now)

integer(i_kind):: lmf           ! number of vertical levels for filtering (generation one)
integer(i_kind):: lmh           ! number of vertical levels for filtering (high generations)
integer(i_kind):: lmf_all       ! number of vertically stacked variables (generation one)
integer(i_kind):: lmh_all       ! number of vertically stacked high generations variabes
integer(i_kind):: kmf           ! number of vertically stacked variables (generation one)
integer(i_kind):: kmh           ! number of vertically stacked high generations variabes


real(r_kind):: lengthx,lengthy,xa0,ya0,xf0,yf0
real(r_kind):: dxf,dyf,dxa,dya

integer(i_kind):: npadx         ! x padding on analysis grid
integer(i_kind):: mpady         ! y padding on analysis grid

integer(i_kind):: ipadx         ! x padding on filter decomposition
integer(i_kind):: jpady         ! y padding on filter deocmposition

!
! Just for standalone test
!
logical:: ldelta
#if 1 
type mg_parameter_type
#endif
contains
procedure,nopass :: init =>init_mg_parameter 
end type mg_parameter_type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                           contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_parameter()
!**********************************************************************!
!                                                                      !
! Initialize ....                                                      !
!                                                                      !
!**********************************************************************!
integer(i_kind):: g

!
! Set number of PEs in x and y directions
!
  namelist /parameters_mgbeta/ mg_ampl01,mg_ampl02,mg_ampl03            &
                              ,mg_weig1,mg_weig2,mg_weig3,mg_weig4      &
                              ,hx,hy,hz,p                               &
                              ,mgbf_line, mgbf_proc                     &
                              ,nxPE,nyPE,im_filt,jm_filt                &
                              ,nm0,mm0                                  &
                              ,lm,lmf,lmh                               &
                              ,ldelta,lquart
!
  open(unit=10,file='mgbeta.nml',status='old',action='read')
  read(10,nml=parameters_mgbeta)
  close(unit=10)
!
!-----------------------------------------------------------------

     nxm = nxPE
     mym = nyPE
!
     im = im_filt
     jm = jm_filt

!-----------------------------------------------------------------
!
!
! For 168 PES
!
!    nxm = 14
!    mym = 12
!
! For 256 PES
!
!
!    nxm =  16
!    mym =  16
!
! For 336 PES
!
!    nxm =  28
!    mym =  12
!
! For 448 PES
!
!    nxm =  28
!    mym =  16
!
!
! For 512 PES
!
!    nxm =  32
!    mym =  16
!
! For 704 PES
!
!    nxm =  32
!    mym =  22
!
! For 768 PES
!
!    nxm =  32
!    mym =  24
!
!
! For 924 PES
!
!    nxm =  28
!    mym =  33
!
! For 1056 PES
!
!    nxm =  32
!    mym =  33
!
! For 1408 PES
!
!    nxm =  32
!    mym =  44
!
! For 1848 PES
!
!    nxm =  56
!    mym =  33
!
! For 2464 PES
!
!    nxm =  56
!    mym =  44


!
! Define maximum number of generations 'gm'
!

      call def_maxgen(nxm,mym,gm)

! Restrict to 4

      if(gm>4) then
        gm=4
      endif
!

!***
!***     Analysis grid
!***

!
! Number of grid intervals on GSI grid for the reduced RTMA domain
! before padding 
!
!clt needed jedi input
  nA_max0 = 1792  
  mA_max0 = 1056  


!
! Number of grid points on the analysis grid after padding
!

!SMALL DOMAIN
!    nm0 = 1792
!    mm0 = 1056
!SMALL DOMAIN

!TEST
!     nm0 = 384
!     mm0 = 384
!TEST
  
    nm = nm0/nxm
    mm = mm0/mym

!***
!***     Filter grid
!***

!    im =  nm
!    jm =  mm

!
! For 168 PES
!
!    im = 120
!    jm =  80

! For 256 PES
!

!    im = 96
!    jm = 64

!    im = 88
!    jm = 56

!
! For 336 PES
!

!    im = 56
!    jm = 80
!
! For 448 PES
!
!    im = 56  
!    jm = 64
!
! For 512 PES
!
!    im = 48
!    jm = 64
!
! For 704 PES
!
!    im = 48  
!    jm = 40
!
! For 768 PES
!
!    im = 48
!    jm = 40
!
! For 924 PES
!
!    im = 56
!    jm = 24
!
! For 1056 PES
!
!    im = 48
!    jm = 24
!
! For 1408 PES
!
!    im = 48
!    jm = 20
!
! For 1848 PES
!
!    im = 28
!    jm = 24
!
! For 2464 PES
!
!    im = 28
!    jm = 20

  im00 = nxm*im
  jm00 = mym*jm

  n0 = 0      ! For now
  m0 = 0      ! For now

  i0 = 1
  j0 = 1

!
! Make sure that nm0 and mm0 and divisibvle with nxm and mym
!
  if(nm*nxm /= nm0 ) then
    write(17,*) 'nm,nxm,nm0=',nm,nxm,nm0
    stop 'nm0 is not divisible by nxm'
  endif
  
  if(mm*mym /= mm0 ) then
    write(17,*) 'mm,mym,mm0=',mm,mym,mm0
    stop 'mm0 is not divisible by mym'
  endif

!
! Set number of processors at higher generations
!

    allocate(ixm(gm))
    allocate(jym(gm))
    allocate(nxy(gm))
    allocate(maxpe_fgen(0:gm))
    allocate(im0(gm))
    allocate(jm0(gm))
    allocate(Fimax(gm))
    allocate(Fjmax(gm))
    allocate(FimaxL(gm))
    allocate(FjmaxL(gm))

    call def_ngens(ixm,gm,nxm)
    call def_ngens(jym,gm,mym)


    do g=1,gm
      nxy(g)=ixm(g)*jym(g)
    enddo

      maxpe_fgen(0)= 0
    do g=1,gm
      maxpe_fgen(g)=maxpe_fgen(g-1)+nxy(g)
    enddo

      maxpe_filt=maxpe_fgen(gm)
      npes_filt=maxpe_filt-nxy(1)

      im0(1)=im00
    do g=2,gm
      im0(g)=im0(g-1)/2
    enddo

      jm0(1)=jm00
    do g=2,gm
      jm0(g)=jm0(g-1)/2
    enddo

    do g=1,gm
      Fimax(g)=im0(g)-im*(ixm(g)-1)
      Fjmax(g)=jm0(g)-jm*(jym(g)-1)
!TEST 
!      write(15,*)'Fimax(',g,')=',Fimax(g)
!      write(15,*)'Fjmax(',g,')=',Fjmax(g)
!TEST 
    enddo

    do g=1,gm
      FimaxL(g)=Fimax(g)/2
      FjmaxL(g)=Fjmax(g)/2
    enddo

!***
!***     Number of variables
!***
!cltmgbf4jedi

  km3 = 6
  km2 = 4

!***
!***     Vertical distribution
!***
 
! lm = 1
!  lm = 50
!  lm05 = lm/2
  km = km3*lm+km2
  lm_all = km3*lm+km2      ! to be deleted

!  lmf = 48
!  lmh = lmf/2
!TEST
!  lmf = lm
!  lmh = lmf
!TEST
  lmf_all = km3*lmf+km2    ! to be deleted
  lmh_all = km3*lmh+km2    ! to be deleted
  kmf = km3*lmf+km2
  kmh = km3*lmh+km2



!***
!*** Filter related parameters
!**
!clt mgbf4jedi to be put into namelist? 
   lengthx = 6.      ! arbitrary chosen scale of the domain
   lengthy = 6.      ! arbitrary chosen scale of the domain


   ib=4
   jb=4

   dxa = lengthx/nm
   dxf = lengthx/im
   nb = 2*dxf/dxa

   dya = lengthy/mm
   dyf = lengthy/jm
   mb = 2*dyf/dya

   xa0 =0.
   ya0 =0.

   xf0=-dxf*0.5
   yf0=-dyf*0.5

  imL=im/2
  jmL=jm/2

!  pasp0=1
!  pasp0 = 5         !  Main
!!  pasp0 = 2.
  pasp01 = mg_ampl01
  pasp02 = mg_ampl02
  pasp03 = mg_ampl03


!TEST  hx=8 
!TEST  hz=8
!TEST  hz=4
!TEST  hz=5
!  hx=6 
!  hy=hx
!  hz=6
!clt mgbf4jedi 
  nh= 6
  nfil = nh + 2

!  p = 4                !  Exponent of Beta function
!  p = 2                !  Exponent of Beta function

  pee2=p*2
    rmom2_1=u1/sqrt(pee2+3)
    rmom2_2=u1/sqrt(pee2+4)
    rmom2_3=u1/sqrt(pee2+5)
    rmom2_4=u1/sqrt(pee2+6)

!----------------------------------------------------------------------
                        end subroutine init_mg_parameter

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_maxgen                          &
!**********************************************************************
!                                                                     !
!  Given number of PEs in x and y direction decides what is the       !
!  maximum number of generations that a multigrid scheme can support  !
!                                                                     !
!                                                    M. Rancic 2020   !
!**********************************************************************
(nxm,mym,gm)
!----------------------------------------------------------------------
implicit none
integer, intent(in):: nxm,mym
integer, intent(out):: gm
integer:: npx,npy,gx,gy

   npx = nxm;  gx=1
   Do 
     npx = (npx + 1)/2
     gx = gx + 1
     if(npx == 1) exit
   end do

   npy = mym;  gy=1
   Do 
     npy = (npy + 1)/2
     gy = gy + 1
     if(npy == 1) exit
   end do

   gm = Min(gx,gy)


!----------------------------------------------------------------------
                        endsubroutine def_maxgen 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_ngens                           &
!*********************************************************************!
!                                                                     !
!  Given number of generations, find number of PEs is s direction     !
!                                                                     !
!                                                    M. Rancic 2020   !
!*********************************************************************!
(nsm,gm,nsm0)
!----------------------------------------------------------------------
implicit none
integer, intent(in):: gm,nsm0
integer, dimension(gm), intent(out):: nsm
integer:: g
!----------------------------------------------------------------------

     nsm(1)=nsm0
   Do g=2,gm
     nsm(g) = (nsm(g-1) + 1)/2
   end do

!----------------------------------------------------------------------
                        endsubroutine def_ngens

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end module mg_parameter
