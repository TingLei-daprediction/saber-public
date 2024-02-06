!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_parameter
!***********************************************************************
!                                                                      !
!       Set resolution, grid and decomposition                         !
!                       - offset version -                             !
!                                                                      !
!  Note:  ixm(1)=nxm, jym(1)=mym                                       !
!                                                                      !
!    If   mod(nxm,2)=0  then  mod(im0,2)=0                             !
!    If   mod(nxm,2)>0  then  mod(im0,8)=0 (for 4 generations)         !
!    (This will keep the right boundary of all decompmisitions at      !
!     same physical location)                                          !
!                                                                      !
!  Modules: kinds, jp_pietc                                            !
!                                                     M. Rancic (2022) !
!***********************************************************************
!clt org use mpi
use kinds, only: i_kind,r_kind
use jp_pietc, only: u1
!use berror, only: mg_ampl0,im_filt,jm_filt
!TEST
!use mpimod, only: nxpe,nype
!TEST

implicit none
type::  mg_parameter_type
!-----------------------------------------------------------------------
!*** 
!*** Namelist parameters
!***
real(r_kind):: mg_ampl01,mg_ampl02,mg_ampl03
real(r_kind):: mg_weig1,mg_weig2,mg_weig3,mg_weig4
integer(i_kind):: mgbf_proc
logical:: mgbf_line
integer(i_kind):: nxPE,nyPE,im_filt,jm_filt
logical:: lquart,lhelm

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
integer(i_kind):: km2_f         ! number of 2d variables for filtering
integer(i_kind):: km3_f         ! number of 3d variables for filtering 
integer(i_kind):: km2_e         ! number of 2d variables for ensemble
integer(i_kind):: km3_e         ! number of 3d variables for ensemble
logical :: l_filt               ! logical flag for filtering or enseble
!integer(i_kind):: km            ! number of vertically stacked all variables (km=km2+lm*km3)
integer(i_kind):: lmf           ! number of vertical levels for filtering (generation one)
integer(i_kind):: lmh           ! number of vertical levels for filtering (high generations)



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
!cltmovedfrom mg_entrymod.f90
integer(i_kind):: km,km2,km3
!cltmoved from type_mg_mppstuff.f90
integer(i_kind):: mype
character(len=5):: c_mype
integer(i_kind):: npes,iTYPE,rTYPE,dTYPE,mpi_comm_comp,ierr,ierror
integer(i_kind):: mpi_comm_work,group_world,group_work

integer(i_kind):: mype_gr,npes_gr

integer(i_kind) my_hgen
integer(i_kind) mype_hgen
logical:: l_hgen
integer(i_kind):: nx,my
!clt moved from *_mg_domain.f90
logical,dimension(2):: Flwest,Fleast,Flnorth,Flsouth
integer(i_kind),dimension(2):: Fitarg_n,Fitarg_e,Fitarg_s,Fitarg_w                         
integer(i_kind),dimension(2):: Fitarg_sw,Fitarg_se,Fitarg_ne,Fitarg_nw

logical,dimension(2):: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne
integer(i_kind),dimension(2):: Fitarg_up 

integer(i_kind):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw


integer(i_kind):: itarg_wA,itarg_eA,itarg_sA,itarg_nA
logical:: lwestA,leastA,lsouthA,lnorthA


integer(i_kind) ix,jy

integer(i_kind),dimension(2):: mype_filt


contains
procedure::init_mg_parameter 
procedure:: init_mg_MPI 
procedure:: finishMPI 
procedure:: barrierMPI 

!clt from mg_domain
procedure:: init_mg_domain
!clt from jp_pbfil.f90
generic :: cholaspect =>cholaspect1,cholaspect2,cholaspect3,cholaspect4
procedure,nopass:: cholaspect1,cholaspect2,cholaspect3,cholaspect4
generic ::  getlinesum=> getlinesum1,getlinesum2,getlinesum3
procedure::   getlinesum1,getlinesum2,getlinesum3
generic :: rbeta=> rbeta1, rbeta2, rbeta3, rbeta4,    &
                   vrbeta1,vrbeta2,vrbeta3,vrbeta4                 
procedure:: rbeta1, rbeta2, rbeta3, rbeta4,    &
                   vrbeta1,vrbeta2,vrbeta3,vrbeta4                 
generic :: rbetaT=>rbeta1t, rbeta2t, rbeta3t, rbeta4t, &
                   vrbeta1t,vrbeta2t,vrbeta3t,vrbeta4t
procedure:: rbeta1t, rbeta2t, rbeta3t, rbeta4t, &
                   vrbeta1t,vrbeta2t,vrbeta3t,vrbeta4t
end type  mg_parameter_type 

interface 
   module subroutine init_mg_MPI(this)
   import :: mg_parameter_type
   class(mg_parameter_type),target :: this
   end subroutine init_mg_MPI
   module subroutine finishMPI(this)
   import :: mg_parameter_type
   class(mg_parameter_type),target :: this
   end subroutine finishMPI
   module subroutine barrierMPI(this)
   import :: mg_parameter_type
   class(mg_parameter_type),target :: this
   end subroutine barrierMPI
   module subroutine init_mg_domain(this) 
     import mg_parameter_type
     class(mg_parameter_type)::this
   end subroutine init_mg_domain
   module subroutine init_domain(this) 
     import mg_parameter_type
     class(mg_parameter_type),target::this
   end subroutine init_domain
   
     module    subroutine init_topology_2d(this)
     import mg_parameter_type
        class(mg_parameter_type),target::this
    end subroutine init_topology_2d
    module     subroutine real_itarg (this,itarg)
     import mg_parameter_type
        class(mg_parameter_type),target::this
        integer(i_kind), intent(inout):: itarg
    end     subroutine real_itarg                     

end interface
!clt from jp_pbfil
!clt from jb_pbfile

interface 

module subroutine cholaspect1(lx,mx, el)                              !  [cholaspect]
use kinds, only: dp=>r_kind
integer,                  intent(in   ):: lx,mx
real(dp),dimension(1,1,lx:mx),intent(inout):: el
!-----------------------------------------------------------------------------
end subroutine cholaspect1
module subroutine cholaspect2(lx,mx, ly,my, el)                        !  [cholaspect]
use kinds, only: dp=>r_kind
integer,                            intent(in   ):: lx,mx, ly,my
real(dp),dimension(2,2,lx:mx,ly:my),intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: tel
end subroutine cholaspect2

module subroutine cholaspect3(lx,mx, ly,my, lz,mz, el)                 !  [cholaspect]
!=============================================================================
! Convert the given field, el, of aspect tensors into the equivalent
! field
! of Cholesky lower-triangular factors of the inverses of the aspect
! tensors.
!=============================================================================
use kinds, only: dp=>r_kind
integer,                                  intent(in   ):: lx,mx, ly,my, lz,mz
real(dp),dimension(3,3,lx:mx,ly:my,lz:mz),intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: tel
end subroutine cholaspect3

module subroutine cholaspect4(lx,mx, ly,my, lz,mz, lw,mw,el)           !  [cholaspect]
use kinds, only: dp=>r_kind
integer,                         intent(in   ):: lx,mx, ly,my, lz,mz, lw,mw
real(dp),dimension(4,4,lx:mx,ly:my,lz:mz,lw:mw),&
                                 intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(4,4):: tel
end subroutine cholaspect4
module subroutine getlinesum1(this,hx,lx,mx, el, ss)                       !  [getlinesum]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                  intent(in   ):: hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
real(dp),dimension(lx:mx),intent(  out):: ss
end subroutine getlinesum1

module subroutine getlinesum2(this,hx,lx,mx, hy,ly,my, el, ss)             !  [getlinesum]
use kinds, only: dp=>r_kind
!=============================================================================
        class(mg_parameter_type)::this
integer,                            intent(in   ):: hx,Lx,mx, &
                                                    hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my),intent(  out):: ss
end subroutine getlinesum2
module subroutine getlinesum3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el, ss)!     [getlinesum]
!=============================================================================
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                  intent(in   ):: hx,Lx,mx, &
                                                          hy,ly,my, &
                                                          hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my,lz:mz),intent(  out):: ss
end subroutine getlinesum3
module subroutine getlinesum4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el, ss)                                                   !  [getlinesum]
use kinds, only: dp=>r_kind
!=============================================================================
        class(mg_parameter_type)::this
integer,                                  intent(in   ):: hx,Lx,mx, &
                                                          hy,ly,my, &
                                                          hz,lz,mz, &
                                                          hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my,lz:mz,Lw:Mw),intent(  out):: ss
end subroutine getlinesum4

module subroutine rbeta1(this,hx,lx,mx, el,ss, a)                               !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(   Lx:Mx),   intent(in   ):: el
real(dp),dimension(   Lx:Mx),   intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
end subroutine rbeta1

module subroutine rbeta2(this,hx,lx,mx, hy,ly,my, el,ss, a)                     !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx, &
                                                            hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
end subroutine rbeta2
module subroutine rbeta3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)            !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                   intent(in   ):: hx,Lx,mx,&
                                                           hy,ly,my,&
                                                           hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz), intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz), intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                              lz-hz:mz+hz),intent(inout):: a
end subroutine rbeta3
module subroutine rbeta4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss,a)!    [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                        intent(in   ):: hx,Lx,mx,&
                                                                hy,ly,my,&
                                                                hz,lz,mz,&
                                                                hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy, &
     lz-hz:mz+hz,lw-hw:mw+hw),                  intent(inout):: a
end subroutine rbeta4
module subroutine rbeta1T(this,hx,lx,mx, el,ss, a)                             !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(  Lx:Mx),    intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
end subroutine rbeta1T
module subroutine rbeta2T(this,hx,lx,mx, hy,ly,my, el,ss, a)                   !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx, &
                                                            hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
end subroutine rbeta2T
module subroutine rbeta3T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)         !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx,&
                                                            hy,ly,my,&
                                                            hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),  intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),  intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                               lz-hz:mz+hz),intent(inout):: a
end subroutine rbeta3T
module subroutine rbeta4T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss, a)                                                   !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                        intent(in   ):: hx,Lx,mx,&
                                                                hy,ly,my,&
                                                                hz,lz,mz,&
                                                                hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                   lz-hz:mz+hz,lw-hw:mw+hw),    intent(inout):: a
end subroutine rbeta4T
module subroutine vrbeta4t(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, &
                                              hw,lw,mw, el,ss, a)!    [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                        intent(in   ):: nv, &
                                                                hx,Lx,mx,&
                                                                hy,ly,my,&
                                                                hz,lz,mz,&
                                                                hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw),                  intent(inout):: a
end subroutine vrbeta4t
module subroutine vrbeta1(this,nv,hx,lx,mx, el,ss, a)                           !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                           intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(1,1, Lx:Mx),    intent(in   ):: el
real(dp),dimension(   Lx:Mx),      intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
end subroutine vrbeta1
module subroutine vrbeta2(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)                 !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx, &
                                                               hy,ly,my
real(dp),dimension(  2,2,Lx:Mx,Ly:My),         intent(in   ):: el
real(dp),dimension(      Lx:Mx,Ly:My),         intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
end subroutine vrbeta2
module subroutine vrbeta3(this,nv, hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)       !  [rbeta]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx,&
                                                               hy,ly,my,&
                                                               hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),     intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),     intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
                                  lz-hz:mz+hz),intent(inout):: a

end subroutine vrbeta3
module subroutine vrbeta4(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss,a)                                                       !  [rbeta]
!=============================================================================
! Vector version of rbeta4 filtering nv fields at once.
!=============================================================================
        class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx,&
                                                               hy,ly,my,&
                                                               hz,lz,mz,&
                                                               hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy, &
     lz-hz:mz+hz,lw-hw:mw+hw),                  intent(inout):: a

end subroutine vrbeta4
module subroutine vrbeta1T(this,nv, hx,lx,mx, el,ss, a)                        !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                           intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),     intent(in   ):: el
real(dp),dimension(   Lx:Mx),      intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
end subroutine vrbeta1T
module subroutine vrbeta2T(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)               !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx, &
                                                               hy,ly,my
real(dp),dimension(   2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(       Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
end subroutine vrbeta2T
module subroutine vrbeta3T(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)     !  [rbetat]
use kinds, only: dp=>r_kind
        class(mg_parameter_type)::this
integer,                                    intent(in   ):: nv,      &
                                                            hx,Lx,mx,&
                                                            hy,ly,my,&
                                                            hz,lz,mz
real(dp),dimension(   3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
real(dp),dimension(       Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
                                lz-hz:mz+hz),intent(inout):: a
end subroutine vrbeta3T
end interface




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                           contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_parameter(this,inputfilename)
!**********************************************************************!
!                                                                      !
! Initialize ....                                                      !
!                                                                      !
!**********************************************************************!
implicit none
class (mg_parameter_type),target:: this
integer(i_kind):: g
character(*):: inputfilename
!*** Namelist parameters as local variable
!***
real(r_kind):: mg_ampl01,mg_ampl02,mg_ampl03
real(r_kind):: mg_weig1,mg_weig2,mg_weig3,mg_weig4
integer(i_kind):: mgbf_proc
logical:: mgbf_line
integer(i_kind):: nxPE,nyPE,im_filt,jm_filt
logical:: lquart,lhelm
logical:: ldelta

integer(i_kind):: lm            ! number of vertical layers
integer(i_kind):: km2_f         ! number of 2d variables for filtering
integer(i_kind):: km3_f         ! number of 3d variables for filtering 
integer(i_kind):: km2_e         ! number of 2d variables for ensemble
integer(i_kind):: km3_e         ! number of 3d variables for ensemble
logical :: l_filt               ! logical flag for filtering or enseble
!integer(i_kind):: km            ! number of vertically stacked all variables (km=km2+lm*km3)
integer(i_kind):: lmf           ! number of vertical levels for filtering (generation one)
integer(i_kind):: lmh           ! number of vertical levels for filtering (high generations)
! Global number of data on Analysis grid
!
integer(i_kind):: nm0        
integer(i_kind):: mm0       


integer(i_kind):: hx,hy,hz
integer(i_kind):: p
!clt include  "type_parameter_locpointer.inc"
!clt include  "type_parameter_point2this.inc"


!
! Set number of PEs in x and y directions
!
  namelist /parameters_mgbeta/ mg_ampl01,mg_ampl02,mg_ampl03            &
                              ,mg_weig1,mg_weig2,mg_weig3,mg_weig4      &
                              ,hx,hy,hz,p                               &
                              ,mgbf_line,mgbf_proc                      &
                              ,lm,lmf,lmh                              &
                              ,km2_f,km3_f,km2_e,km3_e                  &
                              ,l_filt                                   &
                              ,ldelta,lquart,lhelm                      &
                              ,nm0,mm0                                  &
                              ,nxPE,nyPE,im_filt,jm_filt                
!
  open(unit=10,file=inputfilename,status='old',action='read')
  read(10,nml=parameters_mgbeta)
  close(unit=10)
!
!-----------------------------------------------------------------
!clt for safety, copy all namelist loc vars to them of this object
  this%mg_ampl01=mg_ampl01
  this%mg_ampl02=mg_ampl02
  this%mg_ampl03=mg_ampl03            
  this%mg_weig1=mg_weig1
  this%mg_weig2=mg_weig2
  this%mg_weig3=mg_weig3
  this%mg_weig4=mg_weig4
  this%hx=hx
  this%hy=hy
  this%hz=hz
  this%p =p                          
  this%mgbf_line=mgbf_line
  this%mgbf_proc=mgbf_proc          
  this%lm=lm
  this%lmf=lmf
  this%lmh=lmh                    
  this%km2_f=km2_f
  this%km3_f=km3_f
  this%km2_e=km2_e
  this%km3_e=km3_e                
  this%l_filt=l_filt             
  this%ldelta=ldelta
  this%lquart=lquart
  this%lhelm=lhelm 
  this%nm0=nm0
  this%mm0=mm0    
  this%nxPE=nxPE
  this%nyPE=nyPE
  this%im_filt=im_filt
  this%jm_filt=jm_filt                
!clt 


     this%nxm = nxPE
     this%mym = nyPE
!
     this%im = im_filt
     this%jm = jm_filt

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

      call def_maxgen(this%nxm,this%mym,this%gm)

! Restrict to 4

      if(this%gm>4) then
        this%gm=4
      endif
!

!***
!***     Analysis grid
!***

!
! Number of grid intervals on GSI grid for the reduced RTMA domain
! before padding 
!
  this%nA_max0 = 1792  
  this%mA_max0 = 1056  


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
  
    this%nm = this%nm0/this%nxm
    this%mm = this%mm0/this%mym

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

  this%im00 = this%nxm*this%im
  this%jm00 = this%mym*this%jm


  this%n0 = 1 
  this%m0 = 1 

  this%i0 = 1
  this%j0 = 1

!
! Make sure that nm0 and mm0 and divisibvle with nxm and mym
!
  if(this%nm*this%nxm /= this%nm0 ) then
    write(17,*) 'nm,nxm,nm0=',this%nm,this%nxm,this%nm0
    stop 'nm0 is not divisible by nxm'
  endif
  
  if(this%mm*this%mym /= this%mm0 ) then
    write(17,*) 'mm,mym,mm0=',this%mm,this%mym,this%mm0
    stop 'mm0 is not divisible by mym'
  endif

!
! Set number of processors at higher generations
!

    allocate(this%ixm(this%gm))
    allocate(this%jym(this%gm))
    allocate(this%nxy(this%gm))
    allocate(this%maxpe_fgen(0:this%gm))
    allocate(this%im0(this%gm))
    allocate(this%jm0(this%gm))
    allocate(this%Fimax(this%gm))
    allocate(this%Fjmax(this%gm))
    allocate(this%FimaxL(this%gm))
    allocate(this%FjmaxL(this%gm))

    call def_ngens(this%ixm,this%gm,this%nxm)
    call def_ngens(this%jym,this%gm,this%mym)


    do g=1,this%gm
      this%nxy(g)=this%ixm(g)*this%jym(g)
    enddo

      this%maxpe_fgen(0)= 0
    do g=1,this%gm
      this%maxpe_fgen(g)=this%maxpe_fgen(g-1)+this%nxy(g)
    enddo

      this%maxpe_filt=this%maxpe_fgen(this%gm)
      this%npes_filt=this%maxpe_filt-this%nxy(1)

      this%im0(1)=this%im00
    do g=2,this%gm
      this%im0(g)=this%im0(g-1)/2
    enddo

      this%jm0(1)=this%jm00
    do g=2,this%gm
      this%jm0(g)=this%jm0(g-1)/2
    enddo

    do g=1,this%gm
      this%Fimax(g)=this%im0(g)-this%im*(this%ixm(g)-1)
      this%Fjmax(g)=this%jm0(g)-this%jm*(this%jym(g)-1)
!TEST 
!      write(15,*)'Fimax(',g,')=',Fimax(g)
!      write(15,*)'Fjmax(',g,')=',Fjmax(g)
!TEST 
    enddo

    do g=1,this%gm
      this%FimaxL(g)=this%Fimax(g)/2
      this%FjmaxL(g)=this%Fjmax(g)/2
    enddo


!***
!*** Filter related parameters
!**
!D   lengthx = 6.      ! arbitrary chosen scale of the domain
!D   lengthy = 6.      ! arbitrary chosen scale of the domain
   this%lengthx = 1.*this%nm     ! arbitrary chosen scale of the domain
   this%lengthy = 1.*this%mm     ! arbitrary chosen scale of the domain


   this%ib=6  !cltthinkdeb what
   this%jb=6

   this%dxa =this%lengthx/this%nm
   this%dxf = this%lengthx/this%im
   this%nb = 2*this%dxf/this%dxa

   this%dya = this%lengthy/this%mm
   this%dyf = this%lengthy/this%jm
   this%mb = 2*this%dyf/this%dya

!D   xa0 =0.
!D   ya0 =0.
   this%xa0 = this%dxa*0.5
   this%ya0 = this%dya*0.5

!D   xf0=-dxf*0.5
!D   yf0=-dyf*0.5
   this%xf0 = this%dxf*0.5
   this%yf0 = this%dyf*0.5

  this%imL=this%im/2
  this%jmL=this%jm/2

!  pasp0=1
!  pasp0 = 5         !  Main
!!  pasp0 = 2.
  this%pasp01 = mg_ampl01
  this%pasp02 = mg_ampl02
  this%pasp03 = mg_ampl03


!TEST  hx=8 
!TEST  hz=8
!TEST  hz=4
!TEST  hz=5
!  hx=6 
!  hy=hx
!  hz=6

  this%nh= 6
  this%nfil = this%nh + 2

!  p = 4                !  Exponent of Beta function
!  p = 2                !  Exponent of Beta function

    this%pee2=this%p*2
    this%rmom2_1=u1/sqrt(this%pee2+3)
    this%rmom2_2=u1/sqrt(this%pee2+4)
    this%rmom2_3=u1/sqrt(this%pee2+5)
    this%rmom2_4=u1/sqrt(this%pee2+6)
 !clt   call this%init_mg_MPI
 !clt   call this%init_mg_domain

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
