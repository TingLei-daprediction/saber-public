submodule(mg_parameter) jp_pbfil
!$$$  submodule documentation block
!                .      .    .                                       .
! module:   jp_pbfil
!   prgmmr: purser           org: NOAA/EMC            date: 2019-03
!
! abstract:  Codes for the beta filters
!
! module history log:
!   2023-04-19  lei     - object-oriented coding
!   2024-02-20  yokota  - refactoring to apply for GSI
!
! Subroutines Included:
!   cholaspect1 -
!   cholaspect2 -
!   cholaspect3 -
!   cholaspect4 -
!   getlinesum1 -
!   getlinesum2 -
!   getlinesum3 -
!   getlinesum4 -
!   rbeta1 -
!   rbeta2 -
!   rbeta3 -
!   rbeta4 -
!   vrbeta4 -
!   rbeta1T -
!   rbeta2T -
!   rbeta3T -
!   rbeta4T -
!   vrbeta4t -
!   vrbeta1 -
!   vrbeta2 -
!   vrbeta3 -
!   vrbeta1T -
!   vrbeta2T -
!   vrbeta3T -
!
! Functions Included:
!
! remarks:
!   The filters invoke the aspect tensor information encoded by the 
!   Cholesky lower-triangular factors, el, of the INVERSE aspect tensors.
!   The routines, "cholaspect", convert (in place) the field of given
!   aspect tensors A to the equivalent cholesky factors of A^(-1).
!   The routines, "getlinesum" precompute the normalization coefficients
!   for each line (row) of the implied matrix form of the beta filter
!   so that the normalized line sum associated with each point of
!   application becomes unity.
!   This makes the application of each filter significantly faster
!   than having to work out the normalization on the fly.
!   Be sure to have run cholaspect, and then getlinesum, prior to applying
!   the beta filters themselves.
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use mpi
use mgbf_kinds, only: dp=>r_kind
use jp_pietc, only: u1
use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
implicit none
! codex debug/develop for new jim's calibrated function
integer,parameter :: nsres_jim_new=50
real(dp),parameter :: deltabi_jim_new=10.0_dp
real(dp),parameter :: om0_jim_new=1.170561_dp
real(dp),parameter :: om1_jim_new=-4.345418_dp
real(dp),parameter :: om2_jim_new=2.583651_dp
real(dp),parameter :: sres2_jim_new(0:nsres_jim_new) = (/ &
0.9999999_dp,0.8390_dp,0.6805_dp,0.5256_dp,0.3757_dp,0.2328_dp, &
0.1004_dp,-0.01556_dp,-0.1031_dp,-0.1345_dp,-0.02307_dp,0.09016_dp, &
0.03778_dp,-0.01302_dp,-0.04396_dp,-0.03795_dp,0.02998_dp,0.03412_dp, &
0.006966_dp,-0.01677_dp,-0.02450_dp,-0.005024_dp,0.02284_dp,0.01201_dp, &
-0.004674_dp,-0.01482_dp,-0.01056_dp,0.01156_dp,0.01199_dp,0.001285_dp, &
-0.008219_dp,-0.009839_dp,0.001731_dp,0.009906_dp,0.004143_dp,-0.003733_dp, &
-0.007635_dp,-0.003235_dp,0.006907_dp,0.005222_dp,-0.0007235_dp,-0.005267_dp, &
-0.004615_dp,0.003483_dp,0.005193_dp,0.001212_dp,-0.003152_dp,-0.004450_dp, &
-0.00001756_dp,0.004429_dp,0.002340_dp /)
real(dp),parameter :: sres3_jim_new(0:nsres_jim_new) = (/ &
0.9999999_dp,0.8637_dp,0.7058_dp,0.5476_dp,0.3934_dp,0.2467_dp, &
0.1123_dp,-0.001744_dp,-0.08045_dp,-0.08895_dp,0.03107_dp,0.03997_dp, &
0.008750_dp,-0.01796_dp,-0.02029_dp,0.007223_dp,0.01323_dp,0.002605_dp, &
-0.007944_dp,-0.006694_dp,0.004073_dp,0.005567_dp,0.0002230_dp,-0.004496_dp, &
-0.002143_dp,0.002710_dp,0.002604_dp,-0.0006152_dp,-0.002734_dp,-0.0004276_dp, &
0.001852_dp,0.001211_dp,-0.0008620_dp,-0.001619_dp,0.0002304_dp,0.001253_dp, &
0.0004811_dp,-0.0008617_dp,-0.0008482_dp,0.0004500_dp,0.0008209_dp,0.00008178_dp, &
-0.0007499_dp,-0.0003567_dp,0.0004787_dp,0.0005061_dp,-0.0001339_dp,-0.0005853_dp, &
-0.00006824_dp,0.0004252_dp,0.0002775_dp /)
! codex debug/develop for new jim's calibrated function (wbfil variant)
! wbfil replaces zbfil's tabulated s-residual (sres) correction with an exact,
! table-free half-span calibration based on the Euler-Maclaurin moment scheme
! (bcofs). These module-scope variables are filled once by inip_jim_new_wbfil
! and read by bfmoms_jim_new_wbfil / hofnm2_jim_new_wbfil / rcalib1_jim_new_wbfil.
integer,parameter :: nbcof_jim_new_wbfil=9
integer,parameter :: np_jim_new_wbfil=6
real(dp) :: bcofs_jim_new_wbfil(0:nbcof_jim_new_wbfil,0:nbcof_jim_new_wbfil)
real(dp) :: op_jim_new_wbfil,nm2x_jim_new_wbfil
integer  :: p2p3_jim_new_wbfil
logical  :: linit_jim_new_wbfil=.false.

contains

!=============================================================================
module subroutine cholaspect1(lx,mx, el)                        ! [cholaspect]
!=============================================================================
! Convert the given field, el, of aspect tensors into the equivalent
! field
! of Cholesky lower-triangular factors of the inverses of the aspect
! tensors.
!=============================================================================
use jp_pmat, only: inv, l1lm
integer,                  intent(in   ):: lx,mx
real(dp),dimension(1,1,lx:mx),intent(inout):: el
!-----------------------------------------------------------------------------
integer :: ix
!=============================================================================
!$omp parallel do private(ix) schedule(static)
do ix=lx,mx; el(1,1,ix)=u1/sqrt(el(1,1,ix)); enddo
!$omp end parallel do
end subroutine cholaspect1
!=============================================================================
module subroutine cholaspect2(lx,mx, ly,my, el)                 ! [cholaspect]
!=============================================================================
! Convert the given field, el, of aspect tensors into the equivalent
! field
! of Cholesky lower-triangular factors of the inverses of the aspect
! tensors.
!=============================================================================
use jp_pmat, only: inv, l1lm
integer,                            intent(in   ):: lx,mx, ly,my
real(dp),dimension(2,2,lx:mx,ly:my),intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(2,2):: tel
integer                :: ix,iy
!=============================================================================
do iy=ly,my; do ix=lx,mx
   tel=el(:,:,ix,iy); call inv(tel); call l1lm(tel,el(:,:,ix,iy))
enddo;       enddo
end subroutine cholaspect2
!=============================================================================
module subroutine cholaspect3(lx,mx, ly,my, lz,mz, el)          ! [cholaspect]
!=============================================================================
! Convert the given field, el, of aspect tensors into the equivalent
! field
! of Cholesky lower-triangular factors of the inverses of the aspect
! tensors.
!=============================================================================
use jp_pmat, only: inv, l1lm
integer,                                  intent(in   ):: lx,mx, ly,my, lz,mz
real(dp),dimension(3,3,lx:mx,ly:my,lz:mz),intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: tel
integer                :: ix,iy,iz
!=============================================================================
do iz=lz,mz; do iy=ly,my; do ix=lx,mx
   tel=el(:,:,ix,iy,iz); call inv(tel); call l1lm(tel,el(:,:,ix,iy,iz))
enddo;       enddo;       enddo
end subroutine cholaspect3
!=============================================================================
module subroutine cholaspect4(lx,mx, ly,my, lz,mz, lw,mw,el)    ! [cholaspect]
!=============================================================================
! Convert the given field, el, of aspect tensors into the equivalent
! field
! of Cholesky lower-triangular factors of the inverses of the aspect
! tensors.
!=============================================================================
use jp_pmat, only: inv, l1lm
integer,                         intent(in   ):: lx,mx, ly,my, lz,mz, lw,mw
real(dp),dimension(4,4,lx:mx,ly:my,lz:mz,lw:mw),&
                                 intent(inout):: el
!-----------------------------------------------------------------------------
real(dp),dimension(4,4):: tel
integer                :: ix,iy,iz,iw
!=============================================================================
do iw=lw,mw; do iz=lz,mz; do iy=ly,my; do ix=lx,mx
   tel=el(:,:,ix,iy,iz,iw); call inv(tel); call l1lm(tel,el(:,:,ix,iy,iz,iw))
enddo;       enddo;       enddo;       enddo
end subroutine cholaspect4

!=============================================================================
module subroutine getlinesum1(this,hx,lx,mx, el, ss)            ! [getlinesum]
!=============================================================================
! Get inverse of the line-sum of the matrix representing the
! unnormalized
! beta function with aspect tensor pasp=(el*el^T)^(-1), and invert the
! result 
! so it can be used subsequently in the normalized version of this
! filter.
!=============================================================================
class(mg_parameter_type)::this
integer,                  intent(in   ):: hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
real(dp),dimension(lx:mx),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter:: eps=1.e-12
real(dp)          :: s,rr,rrc,exx,x
integer           :: ix,gxl,gxm,gx
!=============================================================================
!clt  write(6,*)'thinkdebss Lx,MX = ',Lx, ' ',Mx
do ix=Lx,Mx
   s=0
   exx=el(1,1,ix)*this%rmom2_1
   x=u1/exx
   gxl=ceiling(-x+eps); gxm=floor( x-eps)
   if(gxl<-hx.or.gxm>hx)&
        stop 'In getlinesum1; filter reach fx becomes too large for hx'
   do gx=gxl,gxm
      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      s=s+rrc**this%p
   enddo
   ss(ix)=u1/s
!clt   write(6,*)'thinkdebss is ',ss(ix)
enddo
end subroutine getlinesum1
module subroutine getlinesum1d(this,hx,lx,mx, el, ss)            ! [getlinesum]
!=============================================================================
!clt from getlinesum1, just reduce e1 to a 1d array
! Get inverse of the line-sum of the matrix representing the
! unnormalized
! beta function with aspect tensor pasp=(el*el^T)^(-1), and invert the
! result 
! so it can be used subsequently in the normalized version of this
! filter.
!=============================================================================
class(mg_parameter_type)::this
integer,                  intent(in   ):: hx,Lx,mx
real(dp),dimension(Lx:Mx),intent(in   ):: el
real(dp),dimension(lx:mx),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter:: eps=1.e-12
real(dp)          :: s,rr,rrc,exx,x
integer           :: ix,gxl,gxm,gx
!=============================================================================
!clt  write(6,*)'thinkdebss Lx,MX = ',Lx, ' ',Mx
do ix=Lx,Mx
   s=0
   exx=el(ix)*this%rmom2_1
   x=u1/exx
   gxl=ceiling(-x+eps); gxm=floor( x-eps)
   if(gxl<-hx.or.gxm>hx) then
        write(error_unit,*) 'thinkdeb7777 exx =',exx,' ',this%rmom2_1,' ',hx,' ',el(ix)
        call flush(error_unit)
        write(error_unit,*) 'In getlinesum1dxx; filter reach fx becomes too large for hx'
        call flush(error_unit)
        stop 'In getlinesum1d; filter reach becomes too large for hy'
   endif
   do gx=gxl,gxm
      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      s=s+rrc**this%p
   enddo
   ss(ix)=u1/s
!clt   write(6,*)'thinkdebss is ',ss(ix)
enddo
end subroutine getlinesum1d
!=============================================================================
module subroutine getlinesum2(this,hx,lx,mx, hy,ly,my, el, ss)  ! [getlinesum]
!=============================================================================
class(mg_parameter_type)::this
integer,                            intent(in   ):: hx,Lx,mx, &
                                                    hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter     :: eps=1.e-12
real(dp),dimension(2,2):: tel
real(dp)               :: s,rr,rrx,rrc,exx,eyy,eyx,x,y,xc
integer                :: ix,gx,gxl,gxm
integer                :: iy,gy,gyl,gym
!=============================================================================
do iy=Ly,My; do ix=Lx,Mx
   s=0
   tel=el(:,:,ix,iy)*this%rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   gyl=ceiling(-y+eps); gym=floor( y-eps)
   if(gyl<-hy.or.gym>hy)&
        stop 'In getlinesum2; filter reach becomes too large for hy'
   do gy=gyl,gym
      y=gy; xc=-y*eyx
      rrx=(y*eyy)**2; x=sqrt(u1-rrx)
      gxl=ceiling((xc-x)/exx+eps); gxm=floor((xc+x)/exx-eps)
      if(gxl<-hx.or.gxm>hx)&
           stop 'In getlinesum2; filter reach becomes too large for hx'
      do gx=gxl,gxm
         x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         s=s+rrc**this%p
      enddo! gx
   enddo! gy
   ss(ix,iy)=u1/s
enddo;  enddo! ix, iy
end subroutine getlinesum2
!=============================================================================
module subroutine getlinesum3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el, ss) ! [getlinesum]
!=============================================================================
class(mg_parameter_type)::this
integer,                                  intent(in   ):: hx,Lx,mx, &
                                                          hy,ly,my, &
                                                          hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my,lz:mz),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter     :: eps=1.e-12
real(dp),dimension(3,3):: tel
real(dp)               :: s,rr,rrx,rry,rrc,&
                          exx,eyy,ezz,eyx,ezx,ezy, x,y,z,xc,yc
integer                :: ix,gx,gxl,gxm
integer                :: iy,gy,gyl,gym
integer                :: iz,gz,gzl,gzm
!=============================================================================
ss=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   s=0
   tel=el(:,:,ix,iy,iz)*this%rmom2_3
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3)
   eyx=tel(2,1); ezx=tel(3,1)
   ezy=tel(3,2)
   z=u1/ezz
   gzl=ceiling(-z+eps); gzm=floor( z-eps)
   if(gzl<-hz.or.gzm>hz)&
        stop 'In getlinesum3; filter reach becomes too large for hz'
   do gz=gzl,gzm
      z=gz;           yc=-z*ezy
      rry=(z*ezz)**2; y =sqrt(u1-rry)
      gyl=ceiling((yc-y)/eyy+eps); gym=floor((yc+y)/eyy-eps)
      if(gyl<-hy.or.gym>hy)&
           stop 'In getlinesum3; filter reach becomes too large for hy'
      do gy=gyl,gym
         y=gy;                  xc=-y*eyx-z*ezx
         rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
         gxl=ceiling((xc-x)/exx+eps); gxm=floor((xc+x)/exx-eps)
         if(gxl<-hx.or.gxm>hx)&
              stop 'In getlinesum3; filter reach becomes too large for hx'
         do gx=gxl,gxm
            x=gx
            rr=rrx+(x*exx-xc)**2; rrc=u1-rr
            s=s+rrc**this%p
         enddo! gx
      enddo! gy
   enddo! gz
   ss(ix,iy,iz)=u1/s
enddo; enddo; enddo! ix, iy, iz
end subroutine getlinesum3
!=============================================================================
module subroutine getlinesum4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el, ss)                                                    ! [getlinesum]
!=============================================================================
class(mg_parameter_type)::this
integer,                                  intent(in   ):: hx,Lx,mx, &
                                                          hy,ly,my, &
                                                          hz,lz,mz, &
                                                          hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    lx:mx,ly:my,lz:mz,Lw:Mw),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter     :: eps=1.e-12
real(dp),dimension(4,4):: tel
real(dp)               :: s,rr,rrx,rry,rrz,rrc, &
                          exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz, x,y,z,w,&
                          xc,yc,zc
integer                :: ix,gx,gxl,gxm
integer                :: iy,gy,gyl,gym
integer                :: iz,gz,gzl,gzm
integer                :: iw,gw,gwl,gwm
!=============================================================================
ss=0
do iw=Lw,Mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   s=0
   tel=el(:,:,ix,iy,iz,iw)*this%rmom2_4
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3); eww=tel(4,4)
   eyx=tel(2,1); ezx=tel(3,1); ewx=tel(4,1)
   ezy=tel(3,2); ewy=tel(4,2)
   ewz=tel(4,3)
   w=u1/eww
   gwl=ceiling(-w+eps); gwm=floor( w-eps)
   if(gwl<-hw.or.gwm>hw)&
        stop 'In getlinesum4; filter reach becomes too large for hw'
   do gw=gwl,gwm
      w=gw;           zc=-w*ewz
      rrz=(w-eww)**2; z =sqrt(u1-rrz)
      gzl=ceiling((zc-z)/ezz+eps); gzm=floor((zc+z)/ezz-eps)
      if(gzl<-hz.or.gzm>hz)&
           stop 'In getlinesum4; filter reach becomes too large for hz'
      do gz=gzl,gzm
         z=gz;                  yc=-z*ezy-w*ewy
         rry=rrz+(z*ezz-zc)**2; y =sqrt(u1-rry)
         gyl=ceiling((yc-y)/eyy+eps); gym=floor((yc+y)/eyy-eps)
         if(gyl<-hy.or.gym>hy)&
              stop 'In getlinesum4; filter reach becomes too large for hy'
         do gy=gyl,gym
            y=gy;                  xc=-y*eyx-z*ezx-w*ewx
            rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
            gxl=ceiling((xc-x)/exx+eps); gxm=floor((xc+x)/exx-eps)
            if(gxl<-hx.or.gxm>hx)&
                 stop 'In getlinesum4; filter reach becomes too large for hx'
            do gx=gxl,gxm
               x=gx
               rr=rrx+(x*exx-xc)**2; rrc=u1-rr
               s=s+rrc**this%p
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
   ss(ix,iy,iz,iw)=u1/s
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
end subroutine getlinesum4

!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rcalib1_jim_new(this,hx,Lx,mx,Lbx,mbx,as,xLb,xmb,el,hxm)
!=============================================================================
class(mg_parameter_type)::this
integer,                      intent(in   ):: hx,Lx,mx
logical,                      intent(in   ):: Lbx,mbx
real(dp),dimension(Lx:Mx),    intent(in   ):: as
real(dp),                     intent(  out):: xLb,xmb
real(dp),dimension(0:1,Lx:Mx),intent(  out):: el
integer,dimension(Lx:Mx),     intent(  out):: hxm
real(dp),dimension(-hx:hx)                :: fs
real(dp)                                  :: b,exx,f,r,rc,rrc,s,x
real(dp)                                  :: rpp3o2_jim_new
integer                                   :: ib,ix,ixp,ixm,gx,gxm,gxn,Lxmix,mxmix
!=============================================================================
xLb=0.0_dp
xmb=0.0_dp
if(Lbx .and. as(Lx)>0.0_dp) xLb=u1/sqrt(as(Lx))
if(mbx .and. as(mx)>0.0_dp) xmb=u1/sqrt(as(mx))
rpp3o2_jim_new=sqrt(real(this%p,dp)+1.5_dp)
do ix=Lx,mx
   b=sqrt(max(as(ix),tiny(1.0_dp)))
   s=rpp3o2_jim_new*b
   r=b*deltabi_jim_new
   ib=int(r)
   if(ib<nsres_jim_new)then
      r=r-ib; rc=u1-r
      if    (this%p==2)then
         s=s+rc*sres2_jim_new(ib)+r*sres2_jim_new(ib+1)
      elseif(this%p==3)then
         s=s+rc*sres3_jim_new(ib)+r*sres3_jim_new(ib+1)
      endif
   endif
   exx=u1/max(s,tiny(1.0_dp))
   el(1,ix)=exx
   gxm=floor(u1/exx)
   hxm(ix)=gxm
   fs(-gxm:gxm)=0.0_dp
   fs(0)=u1
   do gx=-gxm,-1
      rrc=u1-(gx*exx)**2
      f=rrc**this%p
      fs(-gx)=f
      fs(gx)=f
   enddo
   if(Lbx)then
      Lxmix=Lx-ix
      gxn=gxm+Lxmix
      do gx=1,gxn
         x=gx*xLb
         if(x>0.0_dp)x=om0_jim_new+x*(om1_jim_new+x*om2_jim_new)
         ixm=Lxmix-gx
         ixp=Lxmix+gx-1
         fs(ixp)=fs(ixp)+x*fs(ixm)
         fs(ixm)=0.0_dp
      enddo
   endif
   if(mbx)then
      mxmix=mx-ix
      gxn=gxm-mxmix
      do gx=1,gxn
         x=gx*xmb
         if(x>0.0_dp)x=om0_jim_new+x*(om1_jim_new+x*om2_jim_new)
         ixp=mxmix+gx
         ixm=mxmix-gx+1
         fs(ixm)=fs(ixm)+x*fs(ixp)
         fs(ixp)=0.0_dp
      enddo
   endif
   el(0,ix)=u1/sqrt(sum(fs(-gxm:gxm)**2))
enddo
end subroutine rcalib1_jim_new

!=============================================================================
module subroutine rbeta1(this,hx,lx,mx, el,ss, a)                    ! [rbeta]
!=============================================================================
! Perform a radial beta-function filter in 1D.
! It averages the surrounding density values, and so preserves the value
! (in its target region) when presented with a constant-density input
! field.
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx.
! The output data occupy the central region
! Lx <= ix <= Mx.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(   Lx:Mx),   intent(in   ):: el
real(dp),dimension(   Lx:Mx),   intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter             :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: x,tb,s,rr,rrc,frow,exx
integer                        :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   tb=0; s=ss(ix)
   exx=el(ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      tb=tb+frow*a(jx)
   enddo
   b(ix)=tb
enddo
a=b
end subroutine rbeta1
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rbeta1_jim_new(this,hx,lx,mx, el, a)
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: tb,exx,rrc
integer                        :: gx,ix,ixp,ixm
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(ix)
   do gx=ceiling(-u1/exx),-1
      ixp=ix+gx
      ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**this%p*(a(ixp)+a(ixm))
   enddo
   b(ix)=tb*el(0,ix)
enddo
a=b
end subroutine rbeta1_jim_new
module subroutine rbeta3d_1(this,nz,hx,lx,mx, el,ss, a)                    ! [rbeta]
!=============================================================================
!clt modified from rbeta1 to treat files of vertical dimension nz
! Perform a radial beta-function filter in 1D.
! It averages the surrounding density values, and so preserves the value
! (in its target region) when presented with a constant-density input
! field.
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx.
! The output data occupy the central region
! Lx <= ix <= Mx.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: nz,hx,Lx,mx
real(dp),dimension(nz, Lx:Mx),   intent(in   ):: el
real(dp),dimension(nz, Lx:Mx),   intent(in   ):: ss
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter             :: eps=1.e-12
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                       :: x,tb,s,rr,rrc,frow,exx
integer                        :: ix,jx,gx,k
!=============================================================================
b=0
do k=1,nz 
do ix=Lx,Mx
   tb=0; s=ss(k,ix)
   exx=el(k,ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      tb=tb+frow*a(k,jx)
   enddo
   b(k,ix)=tb
enddo
enddo
a=b
end subroutine rbeta3d_1
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rbeta3d_1_jim_new(this,nz,hx,lx,mx, el, a)
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
real(dp),dimension(0:1,nz,Lx:Mx),  intent(in   ):: el
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                          :: tb,exx,rrc
integer                           :: gx,ix,ixp,ixm,k
!=============================================================================
b=0
do k=1,nz
do ix=Lx,Mx
   exx=el(1,k,ix)
   tb=a(k,ix)
   do gx=ceiling(-u1/exx),-1
      ixp=ix+gx
      ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**this%p*(a(k,ixp)+a(k,ixm))
   enddo
   b(k,ix)=tb*el(0,k,ix)
enddo
enddo
a=b
end subroutine rbeta3d_1_jim_new
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rflip1_jim_new(this,hx,lx,mx,Lb,mb,xLb,xmb,a)
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
logical,                        intent(in   ):: Lb,mb
real(dp),                       intent(in   ):: xLb,xmb
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp)                                    :: r
integer                                     :: gx,ixp,ixm,Lxm,mxp
!=============================================================================
if(Lb)then
   Lxm=Lx-1
   do gx=1,hx
      ixm=Lx-gx
      ixp=Lxm+gx
      r=gx*xLb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(ixm)=r*a(ixp)
   enddo
endif
if(mb)then
   mxp=mx+1
   do gx=1,hx
      ixp=mx+gx
      ixm=mxp-gx
      r=gx*xmb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(ixp)=r*a(ixm)
   enddo
endif
end subroutine rflip1_jim_new
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rflip3d_1_jim_new(this,nz,hx,lx,mx,Lb,mb,xLb,xmb,a)
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
logical,                           intent(in   ):: Lb,mb
real(dp),                          intent(in   ):: xLb,xmb
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp)                                       :: r
integer                                        :: gx,ixp,ixm,Lxm,mxp
!=============================================================================
if(Lb)then
   Lxm=Lx-1
   do gx=1,hx
      ixm=Lx-gx
      ixp=Lxm+gx
      r=gx*xLb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(:,ixm)=r*a(:,ixp)
   enddo
endif
if(mb)then
   mxp=mx+1
   do gx=1,hx
      ixp=mx+gx
      ixm=mxp-gx
      r=gx*xmb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(:,ixp)=r*a(:,ixm)
   enddo
endif
end subroutine rflip3d_1_jim_new
!=============================================================================
module subroutine rbeta2(this,hx,lx,mx, hy,ly,my, el,ss, a)          ! [rbeta]
!=============================================================================
! Perform a radial beta-function filter in 2D.
! It averages the surrounding density values, and so preserves the value
! (in its target region) when presented with a constant-density input
! field.
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx, Ly-hy <= Jy <= my+hy
! The output data occupy the central region
! Lx <= ix <= Mx, Ly <= iy <= My.
!=============================================================================
class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx, &
                                                            hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy):: b
real(dp),dimension(2,2)                    :: tel
real(dp)                                   :: tb,s,rr,rrx,rrc,&
                                              frow,exx,eyy,eyx,x,y,xc
integer                                    :: ix,jx,gx
integer                                    :: iy,jy,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy)
   tel=el(:,:,ix,iy)*this%rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy;       y=gy; xc=-y*eyx
      rrx=(y*eyy)**2;       x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**this%p
         tb=tb+frow*a(jx,jy)
      enddo! gx
   enddo! gy
   b(ix,iy)=tb
enddo; enddo! ix, iy
a=b
end subroutine rbeta2
!=============================================================================
module subroutine rbeta3(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a) ! [rbeta]
!=============================================================================
! Perform a radial beta-function filter in 3D.
! It averages the surrounding density values, and so preserves the value
! (in its target region) when presented with a constant-density input
! field.
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx, Ly-hy <= Jy <= my+hy, Lz-hz <= Jz <= mz+hz
! The output data occupy the central region
! Lx <= ix <= Mx, Ly <= iy <= My, Lz <= iz <= Mz.
!=============================================================================
class(mg_parameter_type)::this
integer,                                   intent(in   ):: hx,Lx,mx,&
                                                           hy,ly,my,&
                                                           hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz), intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz), intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                              lz-hz:mz+hz),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                                     :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz):: b
real(dp),dimension(3,3)                                :: tel
real(dp):: s,tb,rr,rrx,rry,rrc,frow,&
           exx,eyy,ezz,eyx,ezx,ezy,x,y,z,xc,yc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy,iz)
   tel=el(:,:,ix,iy,iz)*this%rmom2_3
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3)
   eyx=tel(2,1); ezx=tel(3,1); ezy=tel(3,2)
   z=u1/ezz
   do gz=ceiling(-z+eps),floor( z-eps)
      jz=iz+gz; z=gz; yc=-z*ezy
      rry=(z*ezz)**2; y =sqrt(u1-rry)
      do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
         jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx
         rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
         do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
            jx=ix+gx; x=gx
            rr=rrx+(x*exx-xc)**2; rrc=u1-rr
            frow=s*rrc**this%p
            tb=tb+frow*a(jx,jy,jz)
         enddo! gx
      enddo! gy
   enddo! gz
   b(ix,iy,iz)=tb
enddo;   enddo;    enddo! ix, iy, iz
a=b
end subroutine rbeta3
!=============================================================================
module subroutine rbeta4(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss,a) ! [rbeta]
!=============================================================================
! Perform a radial beta-function filter in 4D.
! It averages the surrounding density values, and so preserves the value
! (in its target region) when presented with a constant-density input
! field.
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx, Ly-hy <= Jy <= my+hy, Lz-hz <= Jz <= mz+hz, 
! Lw-hw <= Jw <= mw+hw
! The output data occupy the central region
! Lx <= ix <= Mx, Ly <= iy <= My, Lz <= iz <= Mz, Lw <= iw <= Mw.
!=============================================================================
class(mg_parameter_type)::this
integer,                                        intent(in   ):: hx,Lx,mx,&
                                                                hy,ly,my,&
                                                                hz,lz,mz,&
                                                                hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy, &
     lz-hz:mz+hz,lw-hw:mw+hw),                  intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw)              :: b
real(dp),dimension(4,4)                    :: tel
real(dp):: s,tb,rr,rrx,rry,rrz,rrc,frow,&
           exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,x,y,z,w,xc,yc,zc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
integer :: iw,jw,gw
!=============================================================================
b=0
do iw=lw,mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy,iz,iw)
   tel=el(:,:,ix,iy,iz,iw)*this%rmom2_4
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3); eww=tel(4,4)
   eyx=tel(2,1); ezx=tel(3,1); ewx=tel(4,1)
   ezy=tel(3,2); ewy=tel(4,2)
   ewz=tel(4,3)
   w=u1/eww
   do gw=ceiling(-w+eps),floor( w-eps)
      jw=iw+gw; w=gw; zc=-w*ewz
      rrz=(w*eww)**2; z =sqrt(u1-rrz)
      do gz=ceiling((zc-z)/ezz+eps),floor((zc+z)/ezz-eps)
         jz=iz+gz; z=gz;        yc=-z*ezy-w*ewy
         rry=rrz+(z*ezz-zc)**2; y =sqrt(u1-rry)
         do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
            jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx-w*ewx
            rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
            do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
               jx=ix+gx; x=gx
               rr=rrx+(x*exx-xc)**2; rrc=u1-rr
               frow=s*rrc**this%p
               tb=tb+frow*a(jx,jy,jz,jw)
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
   b(ix,iy,iz,iw)=tb
enddo;   enddo;   enddo;   enddo! ix, iy, iz, iw
a=b
end subroutine rbeta4

!=============================================================================
! Vector versions of the above routines:
!=============================================================================
module subroutine vrbeta4(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss,a)                                                        ! [rbeta]
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
!-----------------------------------------------------------------------------
real(dp),parameter                                     :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw)                          :: b
real(dp),dimension(nv)                                 :: tb
real(dp),dimension(4,4)                                :: tel
real(dp):: s,rr,rrx,rry,rrz,rrc,frow,&
           exx,eyy,ezz,eww, eyx,ezx,ewx, ezy,ewy, ewz,&
           x,y,z,w,xc,yc,zc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
integer :: iw,jw,gw
!=============================================================================
b=0
do iw=lw,mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy,iz,iw)
   tel=el(:,:,ix,iy,iz,iw)*this%rmom2_4
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3); eww=tel(4,4)
   eyx=tel(2,1); ezx=tel(3,1); ewx=tel(4,1)
   ezy=tel(3,2); ewy=tel(4,2)
   ewz=tel(4,3)
   w=u1/eww
   do gw=ceiling(-w+eps),floor( w-eps)
      jw=iw+gw; w=gw; zc=-w*ewz
      rrz=(w*eww)**2; z =sqrt(u1-rrz)
      do gz=ceiling((zc-z)/ezz+eps),floor((zc+z)/ezz-eps)
         jz=iz+gz; z=gz;        yc=-z*ezy-w*ewy
         rry=rrz+(z*ezz-zc)**2; y =sqrt(u1-rry)
         do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
            jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx-w*ewx
            rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
            do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
               jx=ix+gx; x=gx
               rr=rrx+(x*exx-xc)**2; rrc=u1-rr
               frow=s*rrc**this%p
               tb=tb+frow*a(:,jx,jy,jz,jw)
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
   b(:,ix,iy,iz,iw)=tb
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
a=b
end subroutine vrbeta4

!=============================================================================
module subroutine rbeta1T(this,hx,lx,mx, el,ss, a)                  ! [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 1D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= mx+hx.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(  Lx:Mx),    intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter             :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: ta,s,rr,rrc,frow,exx,x
integer                        :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   ta=a(ix); s=ss(ix)
   exx=el(1,1,ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      b(jx)=b(jx)+frow*ta
   enddo
enddo
a=b
end subroutine rbeta1t
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rbeta1T_jim_new(this,hx,lx,mx, el, a)
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: ta,exx,rrc,tafrow
integer                        :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   ta=a(ix)*el(0,ix)
   exx=el(1,ix)
   b(ix)=b(ix)+ta
   do gx=ceiling(-u1/exx),-1
      jx=ix+gx
      rrc=u1-(gx*exx)**2
      tafrow=ta*rrc**this%p
      b(jx)=b(jx)+tafrow
      b(ix-gx)=b(ix-gx)+tafrow
   enddo
enddo
a=b
end subroutine rbeta1T_jim_new
module subroutine rbeta3d_1T(this,nz,hx,lx,mx, el,ss, a)                  ! [rbetat]
!clt modified from rbeta1T to add a vertical dimension
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 1D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= mx+hx.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   )::nz, hx,Lx,mx
real(dp),dimension(nz,Lx:Mx),  intent(in   ):: el
real(dp),dimension(nz,  Lx:Mx),    intent(in   ):: ss
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter             :: eps=1.e-12
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                       :: ta,s,rr,rrc,frow,exx,x
integer                        :: ix,jx,gx,k
!=============================================================================
b=0
do k=1,nz
do ix=Lx,Mx
   ta=a(k,ix); s=ss(k,ix)
   exx=el(k,ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      b(k,jx)=b(k,jx)+frow*ta
   enddo
enddo
enddo
a=b
end subroutine rbeta3d_1t
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rbeta3d_1T_jim_new(this,nz,hx,lx,mx, el, a)
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
real(dp),dimension(0:1,nz,Lx:Mx),  intent(in   ):: el
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                          :: ta,exx,rrc,tafrow
integer                           :: ix,jx,gx,k
!=============================================================================
b=0
do k=1,nz
do ix=Lx,Mx
   ta=a(k,ix)*el(0,k,ix)
   exx=el(1,k,ix)
   b(k,ix)=b(k,ix)+ta
   do gx=ceiling(-u1/exx),-1
      jx=ix+gx
      rrc=u1-(gx*exx)**2
      tafrow=ta*rrc**this%p
      b(k,jx)=b(k,jx)+tafrow
      b(k,ix-gx)=b(k,ix-gx)+tafrow
   enddo
enddo
enddo
a=b
end subroutine rbeta3d_1T_jim_new
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rflip1T_jim_new(this,hx,lx,mx,Lb,mb,xLb,xmb,a)
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
logical,                        intent(in   ):: Lb,mb
real(dp),                       intent(in   ):: xLb,xmb
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp)                                    :: r
integer                                     :: gx,ixp,ixm,Lxm,mxp
!=============================================================================
if(Lb)then
   Lxm=Lx-1
   do gx=1,hx
      ixm=Lx-gx
      ixp=Lxm+gx
      r=gx*xLb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(ixp)=a(ixp)+r*a(ixm)
      a(ixm)=0.0_dp
   enddo
endif
if(mb)then
   mxp=mx+1
   do gx=1,hx
      ixp=mx+gx
      ixm=mxp-gx
      r=gx*xmb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(ixm)=a(ixm)+r*a(ixp)
      a(ixp)=0.0_dp
   enddo
endif
end subroutine rflip1T_jim_new
!=============================================================================
! codex debug/develop for new jim's calibrated function
module subroutine rflip3d_1T_jim_new(this,nz,hx,lx,mx,Lb,mb,xLb,xmb,a)
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
logical,                           intent(in   ):: Lb,mb
real(dp),                          intent(in   ):: xLb,xmb
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp)                                       :: r
integer                                        :: gx,ixp,ixm,Lxm,mxp
!=============================================================================
if(Lb)then
   Lxm=Lx-1
   do gx=1,hx
      ixm=Lx-gx
      ixp=Lxm+gx
      r=gx*xLb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(:,ixp)=a(:,ixp)+r*a(:,ixm)
      a(:,ixm)=0.0_dp
   enddo
endif
if(mb)then
   mxp=mx+1
   do gx=1,hx
      ixp=mx+gx
      ixm=mxp-gx
      r=gx*xmb
      r=om0_jim_new+r*(om1_jim_new+r*om2_jim_new)
      a(:,ixm)=a(:,ixm)+r*a(:,ixp)
      a(:,ixp)=0.0_dp
   enddo
endif
end subroutine rflip3d_1T_jim_new
!=============================================================================
! codex debug/develop for new jim's calibrated function (wbfil variant)
! Ported from dr-jim-new/wbfil.f90. Same el(0:1,..) coefficient format as the
! zbfil-derived *_jim_new routines, but the half-span is calibrated exactly via
! the moment scheme (inip/bfmoms/nm2ofh/hofnm2) instead of the sres lookup
! tables, and rcalib carries no boundary (flip) treatment (wbfil design).
!=============================================================================
module subroutine inip_jim_new_wbfil(this)
!=============================================================================
! One-time initialization of the wbfil moment-calibration constants from this%p.
! Fills bcofs (Euler-Maclaurin end-correction coefficients), op=1/p, the
! normalized 2nd moment nm2x at half-span 2, and p2p3=2p+3.
!=============================================================================
class(mg_parameter_type)::this
integer :: p,i
integer,dimension(0:nbcof_jim_new_wbfil,0:nbcof_jim_new_wbfil):: inums
integer,dimension(0:nbcof_jim_new_wbfil):: idens
real(dp),parameter :: u2=2.0_dp,u3=3.0_dp,u4=4.0_dp
real(dp) :: u4o3
data inums/2,9*0, 1,2,8*0, -1,10,6,7*0, 1,-7,21,6,6*0, -3,20,-42,60,10,5*0,&
     5,-33,66,-66,55,6,4*0, -691,4550,-9009,8580,-5005,2730,210,3*0,       &
     105,-691,1365,-1287,715,-273,105,6,0,0,                               &
     -3617,23800,-46988,44200,-24310,8840,-2380,680,30,0,                  &
     219335,-1443183,2848860,-2678316,1469650,-529074,135660,-27132,5985,210/
data idens/1,3,15,21,45,33,1365,45,255,1995/
!=============================================================================
p=this%p
u4o3=u4/u3
op_jim_new_wbfil=u1/p
nm2x_jim_new_wbfil=u2/(u4o3**p+u2)! normalized 2nd moment of a half-span-2 filter
p2p3_jim_new_wbfil=p*2+3
bcofs_jim_new_wbfil=inums
do i=0,nbcof_jim_new_wbfil; bcofs_jim_new_wbfil(:,i)=bcofs_jim_new_wbfil(:,i)/idens(i); enddo
linit_jim_new_wbfil=.true.
end subroutine inip_jim_new_wbfil
!=============================================================================
module subroutine bfmoms_jim_new_wbfil(this,h,mom0,mom2,dmom0,dmom2)
!=============================================================================
! Exact 0th and 2nd moments (and their h-derivatives) of a beta line filter of
! exponent this%p and half-span h, via the residual-free Euler-Maclaurin scheme.
!=============================================================================
class(mg_parameter_type)::this
real(dp),intent(in ):: h
real(dp),intent(out):: mom0,mom2,dmom0,dmom2
real(dp),dimension(0:np_jim_new_wbfil):: pchoose
real(dp),parameter :: o2=0.5_dp
real(dp) :: ho2,hh,c,dc,q0,q2,enn
integer  :: j,jm,jp,k,n,p
!=============================================================================
p=this%p
pchoose(0)=u1
do j=1,p; jm=j-1; pchoose(j)=(pchoose(jm)*(p-jm))/j; enddo
ho2=h*o2; hh=h*h; n=h; enn=n*n
mom0=0; dmom0=0; mom2=0; dmom2=0
do j=0,p; jp=j+1
   c=pchoose(j)/(-hh)**j; dc=-j*c/ho2
   q0=enn**j;             q2=q0*enn
   do k=0,j;  q0=q0+n*bcofs_jim_new_wbfil(k,j )*enn**k; enddo
   do k=0,jp; q2=q2+n*bcofs_jim_new_wbfil(k,jp)*enn**k; enddo
   mom0  =mom0 +c*q0;    mom2= mom2+ c*q2
   dmom0=dmom0+dc*q0;   dmom2=dmom2+dc*q2
enddo
end subroutine bfmoms_jim_new_wbfil
!=============================================================================
module subroutine nm2ofh_jim_new_wbfil(this,h,nm2,dnm2)
!=============================================================================
! Normalized 2nd moment nm2 (and derivative dnm2) of the beta line filter of
! half-span h.
!=============================================================================
class(mg_parameter_type)::this
real(dp),intent(in ):: h
real(dp),intent(out):: nm2,dnm2
real(dp):: mom0,mom2,dmom0,dmom2,omom0
!=============================================================================
call this%bfmoms_jim_new_wbfil(h,mom0,mom2,dmom0,dmom2)
omom0=u1/mom0
nm2=mom2*omom0
dnm2=(dmom2*mom0-dmom0*mom2)*omom0**2
end subroutine nm2ofh_jim_new_wbfil
!=============================================================================
module subroutine hofnm2_jim_new_wbfil(this,nm2t,h)
!=============================================================================
! Invert nm2ofh: find the half-span h whose single beta filter has normalized
! 2nd moment nm2t. Closed form for h<=2, else Newton iteration.
!=============================================================================
class(mg_parameter_type)::this
real(dp),intent(in ):: nm2t
real(dp),intent(out):: h
integer,parameter :: nit=40
real(dp),parameter:: eps=1.e-12_dp,u2=2.0_dp
real(dp):: nm2,dnm2,r
integer :: it
!=============================================================================
if(nm2t<nm2x_jim_new_wbfil)then
   h=u1/sqrt(u1-(nm2t/(u2*(u1-nm2t)))**op_jim_new_wbfil)
   return
endif
h=sqrt(u1+p2p3_jim_new_wbfil*nm2t)
do it=1,nit
   call this%nm2ofh_jim_new_wbfil(h,nm2,dnm2)
   r=(nm2-nm2t)/dnm2
   h=h-r
   if(abs(r)<eps)return
enddo
end subroutine hofnm2_jim_new_wbfil
!=============================================================================
module subroutine rcalib1_jim_new_wbfil(this,hx,Lx,mx,as,el,hxm)
!=============================================================================
! wbfil 1D calibration: for each aspect as, find the exact half-span whose
! single beta filter has normalized 2nd moment as/2, store its reciprocal in
! el(1,:) and the amplitude normalization in el(0,:). No boundary treatment.
!=============================================================================
class(mg_parameter_type)::this
integer,                      intent(in ):: hx,Lx,mx
real(dp),dimension(    Lx:mx),intent(in ):: as
real(dp),dimension(0:1,Lx:mx),intent(out):: el
integer,dimension(Lx:mx),     intent(out):: hxm
real(dp),dimension(-hx:hx):: fs
real(dp),parameter        :: o2=0.5_dp
real(dp)                  :: exx,f,rrc,s
integer                   :: ix,gx,gxm
!=============================================================================
do ix=Lx,mx
   call this%hofnm2_jim_new_wbfil(as(ix)*o2,s)
   exx=u1/s
   el(1,ix)=exx
   gxm=floor(u1/exx); hxm(ix)=gxm
   fs(-gxm:gxm)=0
   fs(0)=u1
   do gx=-gxm,-1
      rrc=u1-(gx*exx)**2
      f=rrc**this%p; fs(-gx)=f; fs(gx)=f
   enddo
   el(0,ix)=u1/sqrt(sum(fs(-gxm:gxm)**2))
enddo
end subroutine rcalib1_jim_new_wbfil
!=============================================================================
module subroutine rbeta3d_1_jim_new_wbfil(this,nz,hx,lx,mx, el, a)
!=============================================================================
! Direct beta line filter over nz levels (wbfil coefficients). Body matches
! rbeta3d_1_jim_new; only the el coefficients (from rcalib1_jim_new_wbfil) differ.
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
real(dp),dimension(0:1,nz,Lx:Mx),  intent(in   ):: el
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                          :: tb,exx,rrc
integer                           :: gx,ix,ixp,ixm,k
!=============================================================================
b=0
do k=1,nz
do ix=Lx,Mx
   exx=el(1,k,ix)
   tb=a(k,ix)
   do gx=ceiling(-u1/exx),-1
      ixp=ix+gx
      ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**this%p*(a(k,ixp)+a(k,ixm))
   enddo
   b(k,ix)=tb*el(0,k,ix)
enddo
enddo
a=b
end subroutine rbeta3d_1_jim_new_wbfil
!=============================================================================
module subroutine rbeta3d_1T_jim_new_wbfil(this,nz,hx,lx,mx, el, a)
!=============================================================================
! Adjoint of rbeta3d_1_jim_new_wbfil. Body matches rbeta3d_1T_jim_new.
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nz,hx,Lx,mx
real(dp),dimension(0:1,nz,Lx:Mx),  intent(in   ):: el
real(dp),dimension(nz,lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(nz,lx-hx:mx+hx):: b
real(dp)                          :: ta,exx,rrc,tafrow
integer                           :: ix,jx,gx,k
!=============================================================================
b=0
do k=1,nz
do ix=Lx,Mx
   ta=a(k,ix)*el(0,k,ix)
   exx=el(1,k,ix)
   b(k,ix)=b(k,ix)+ta
   do gx=ceiling(-u1/exx),-1
      jx=ix+gx
      rrc=u1-(gx*exx)**2
      tafrow=ta*rrc**this%p
      b(k,jx)=b(k,jx)+tafrow
      b(k,ix-gx)=b(k,ix-gx)+tafrow
   enddo
enddo
enddo
a=b
end subroutine rbeta3d_1T_jim_new_wbfil
!=============================================================================
module subroutine rbeta1_jim_new_wbfil(this,hx,lx,mx, el, a)
!=============================================================================
! wbfil variant of rbeta1_jim_new (1D direct beta line filter). Body is
! identical to rbeta1_jim_new; it is paired with the rcalib1_jim_new_wbfil
! (hofnm2) calibrated coefficients for the jim_new_wbfil path.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: tb,exx,rrc
integer                        :: gx,ix,ixp,ixm
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(ix)
   do gx=ceiling(-u1/exx),-1
      ixp=ix+gx
      ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**this%p*(a(ixp)+a(ixm))
   enddo
   b(ix)=tb*el(0,ix)
enddo
a=b
end subroutine rbeta1_jim_new_wbfil
!=============================================================================
module subroutine rbeta1T_jim_new_wbfil(this,hx,lx,mx, el, a)
!=============================================================================
! Adjoint of rbeta1_jim_new_wbfil. Body matches rbeta1T_jim_new.
!=============================================================================
class(mg_parameter_type)::this
integer,                        intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:Mx),  intent(in   ):: el
real(dp),dimension(lx-hx:mx+hx),intent(inout):: a
real(dp),dimension(lx-hx:mx+hx):: b
real(dp)                       :: ta,exx,rrc,tafrow
integer                        :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   ta=a(ix)*el(0,ix)
   exx=el(1,ix)
   b(ix)=b(ix)+ta
   do gx=ceiling(-u1/exx),-1
      jx=ix+gx
      rrc=u1-(gx*exx)**2
      tafrow=ta*rrc**this%p
      b(jx)=b(jx)+tafrow
      b(ix-gx)=b(ix-gx)+tafrow
   enddo
enddo
a=b
end subroutine rbeta1T_jim_new_wbfil
!=============================================================================
module subroutine rbeta2T(this,hx,lx,mx, hy,ly,my, el,ss, a)        ! [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 2D.
! It conserved "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx, Ly <= iy <= My.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= mx+hx, Ly-hy <= Jy <= my+hy
!=============================================================================
class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx, &
                                                            hy,ly,my
real(dp),dimension(2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy):: b
real(dp),dimension(2,2)                    :: tel
real(dp)                                   :: ta,s,rr,rrx,rrc, &
                                              frow,exx,eyy,eyx,x,y,xc
integer                                    :: ix,jx,gx
integer                                    :: iy,jy,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   ta=a(ix,iy); s=ss(ix,iy)
   tel=el(:,:,ix,iy)*this%rmom2_2 ! sThis el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy; y=gy; xc=-y*eyx
      rrx=(y*eyy)**2; x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**this%p
         b(jx,jy)=b(jx,jy)+frow*ta
      enddo! gx 
   enddo! gy
enddo;  enddo! ix, iy
a=b
end subroutine rbeta2t
!=============================================================================
module subroutine rbeta3T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a) ! [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 3D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx, Ly <= iy <= My, Lz <= iz <= Mz.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= Mx+hx, Ly-hy <= Jy <= My+hy, Lz-hz <= Jz <= Mz+hz.
!=============================================================================
class(mg_parameter_type)::this
integer,                                    intent(in   ):: hx,Lx,mx,&
                                                            hy,ly,my,&
                                                            hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),  intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),  intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                               lz-hz:mz+hz),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                                     :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz):: b
real(dp),dimension(3,3)                                :: tel
real(dp):: ta,s,rr,rrx,rry,rrc,frow,&
           exx,eyy,ezz,eyx,ezx,ezy,x,y,z,xc,yc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   ta=a(ix,iy,iz); s=ss(ix,iy,iz)
   tel=el(:,:,ix,iy,iz)*this%rmom2_3
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3)
   eyx=tel(2,1); ezx=tel(3,1); ezy=tel(3,2)
   z=u1/ezz
   do gz=ceiling(-z+eps),floor( z-eps)
      jz=iz+gz; z=gz; yc=-z*ezy
      rry=(z*ezz)**2; y =sqrt(u1-rry)
      do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
         jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx
         rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
         do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
            jx=ix+gx; x=gx
            rr=rrx+(x*exx-xc)**2; rrc=u1-rr
            frow=s*rrc**this%p
            b(jx,jy,jz)=b(jx,jy,jz)+frow*ta
         enddo! gx
      enddo! gy
   enddo ! gz
enddo;  enddo;  enddo ! ix, iy, iz
a=b
end subroutine rbeta3t
!=============================================================================
module subroutine rbeta4T(this,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss, a)                                                      ! [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 4D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx, Ly <= iy <= My, Lz <= iz <= Mz, Lw <= iw <= Mw.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= Mx+hx, Ly-hy <= Jy <= My+hy, Lz-hz <= Jz <= Mz+hz, 
!     Lw-hw <= Jw <= Mw+hw.
!=============================================================================
class(mg_parameter_type)::this
integer,                                        intent(in   ):: hx,Lx,mx,&
                                                                hy,ly,my,&
                                                                hz,lz,mz,&
                                                                hw,lw,mw
real(dp),dimension(4,4,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: ss
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
                   lz-hz:mz+hz,lw-hw:mw+hw),    intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw)              :: b
real(dp),dimension(4,4)                    :: tel
real(dp):: ta,s,rr,rrx,rry,rrz,rrc,frow,&
           exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,x,y,z,w,xc,yc,zc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
integer :: iw,jw,gw
!=============================================================================
b=0
do iw=Lw,Mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   ta=a(ix,iy,iz,iw); s=ss(ix,iy,iz,iw)
   tel=el(:,:,ix,iy,iz,iw)*this%rmom2_4
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3); eww=tel(4,4)
   eyx=tel(2,1); ezx=tel(3,1); ewx=tel(4,1)
   ezy=tel(3,2); ewy=tel(4,2)
   ewz=tel(4,3)
   z=u1/ezz
   do gw=ceiling(-w+eps),floor( w-eps)
      jw=iw+gw; w=gw; zc=-w*ewz
      rrz=(w*eww)**2; z =sqrt(u1-rrz)
      do gz=ceiling((zc-z)/ezz+eps),floor((zc+z)/ezz-eps)
         jz=iz+gz; z=gz;        yc=-z*ezy-w*ewy
         rry=rrz+(z*ezz-zc)**2; y =sqrt(u1-rry)
         do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
            jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx-w*ewx
            rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
            do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
               jx=ix+gx; x=gx
               rr=rrx+(x*exx-xc)**2; rrc=u1-rr
               frow=s*rrc**this%p
               b(jx,jy,jz,jw)=b(jx,jy,jz,jw)+frow*ta
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
a=b
end subroutine rbeta4t


!=============================================================================
module subroutine vrbeta4t(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, &
                                              hw,lw,mw, el,ss, a)   ! [rbetat]
!=============================================================================
! Vector version of rbeta4t filtering nv fields at once.
!=============================================================================
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
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw)              :: b
real(dp),dimension(nv)                     :: ta
real(dp),dimension(4,4)                    :: tel
real(dp):: s,rr,rrx,rry,rrz,rrc,frow,&
           exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,x,y,z,w,xc,yc,zc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
integer :: iw,jw,gw
!=============================================================================
b=0
do iw=Lw,Mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   ta=a(:,ix,iy,iz,iw); s=ss(ix,iy,iz,iw)
   tel=el(:,:,ix,iy,iz,iw)*this%rmom2_4
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3); eww=tel(4,4)
   eyx=tel(2,1); ezx=tel(3,1); ewx=tel(4,1)
   ezy=tel(3,2); ewy=tel(4,2)
   ewz=tel(4,3)
   z=u1/ezz
   do gw=ceiling(-w+eps),floor( w-eps)
      jw=iw+gw; w=gw; zc=-w*ewz
      rrz=(w*eww)**2; z =sqrt(u1-rrz)
      do gz=ceiling((zc-z)/ezz+eps),floor((zc+z)/ezz-eps)
         jz=iz+gz; z=gz;        yc=-z*ezy-w*ewy
         rry=rrz+(z*ezz-zc)**2; y =sqrt(u1-rry)
         do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
            jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx-w*ewx
            rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
            do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
               jx=ix+gx; x=gx
               rr=rrx+(x*exx-xc)**2; rrc=u1-rr
               frow=s*rrc**this%p
               b(:,jx,jy,jz,jw)=b(:,jx,jy,jz,jw)+frow*ta
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
enddo; enddo; enddo; enddo! ix, iy, iz, iw
a=b
end subroutine vrbeta4t

! Vector versions of the above routines:
!=============================================================================
module subroutine vrbeta1(this,nv,hx,lx,mx, el,ss, a)                ! [rbeta]
!=============================================================================
! Vector version of rbeta1 filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(1,1, Lx:Mx),    intent(in   ):: el
real(dp),dimension(   Lx:Mx),      intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx):: b
real(dp),dimension(nv)            :: tb
real(dp)                          :: x,s,rr,rrc,frow,exx
integer                           :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   tb=0; s=ss(ix)
   exx=el(1,1,ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      tb=tb+frow*a(:,jx)
   enddo
   b(:,ix)=tb
enddo
a=b
end subroutine vrbeta1

!=============================================================================
module subroutine vrbeta2(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)      ! [rbeta]
!=============================================================================
! Vector version of rbeta2 filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx, &
                                                               hy,ly,my
real(dp),dimension(  2,2,Lx:Mx,Ly:My),         intent(in   ):: el
real(dp),dimension(      Lx:Mx,Ly:My),         intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                            :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy):: b
real(dp),dimension(nv)                        :: tb
real(dp),dimension(2,2)                    :: tel
real(dp)                                   :: s,rr,rrx,rrc,&
                                              frow,exx,eyy,eyx,x,y,xc
integer                                    :: ix,jx,gx
integer                                    :: iy,jy,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy)
   tel=el(:,:,ix,iy)*this%rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy;       y=gy; xc=-y*eyx
      rrx=(y*eyy)**2;       x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**this%p
         tb=tb+frow*a(:,jx,jy)
      enddo! gx
   enddo! gy
   b(:,ix,iy)=tb
enddo;   enddo! ix, iy
a=b
end subroutine vrbeta2

!=============================================================================
module subroutine vrbeta3(this,nv, hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a) ! [rbeta]
!=============================================================================
! Vector version of rbeta3 filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx,&
                                                               hy,ly,my,&
                                                               hz,lz,mz
real(dp),dimension(3,3,Lx:Mx,Ly:My,Lz:Mz),     intent(in   ):: el
real(dp),dimension(    Lx:Mx,Ly:My,Lz:Mz),     intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
                                  lz-hz:mz+hz),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                                        :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz):: b
real(dp),dimension(nv)                                    :: tb
real(dp),dimension(3,3)                                   :: tel
real(dp):: s,rr,rrx,rry,rrc,frow,&
           exx,eyy,ezz,eyx,ezx,ezy,x,y,z,xc,yc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tb=0; s=ss(ix,iy,iz)
   tel=el(:,:,ix,iy,iz)*this%rmom2_3
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3)
   eyx=tel(2,1); ezx=tel(3,1); ezy=tel(3,2)
   z=u1/ezz
   do gz=ceiling(-z+eps),floor( z-eps)
      jz=iz+gz; z=gz; yc=-z*ezy
      rry=(z*ezz)**2; y =sqrt(u1-rry)
      do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
         jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx
         rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
         do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
            jx=ix+gx; x=gx
            rr=rrx+(x*exx-xc)**2; rrc=u1-rr
            frow=s*rrc**this%p
            tb=tb+frow*a(:,jx,jy,jz)
         enddo! gx
      enddo! gy
   enddo! gz
   b(:,ix,iy,iz)=tb
enddo;   enddo;    enddo! ix, iy, iz
a=b
end subroutine vrbeta3

! Vector versions of the above routines:
!=============================================================================
module subroutine vrbeta1T(this,nv, hx,lx,mx, el,ss, a)             ! [rbetat]
!=============================================================================
! Vector version of rbeta1t filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                           intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),     intent(in   ):: el
real(dp),dimension(   Lx:Mx),      intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx):: b
real(dp),dimension(nv)            :: ta
real(dp)                          :: s,rr,rrc,frow,exx,x
integer                           :: ix,jx,gx
!=============================================================================
b=0
do ix=Lx,Mx
   ta=a(:,ix); s=ss(ix)
   exx=el(1,1,ix)*this%rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**this%p
      b(:,jx)=b(:,jx)+frow*ta
   enddo
enddo
a=b
end subroutine vrbeta1t
!=============================================================================
module subroutine vrbeta2T(this,nv,hx,lx,mx, hy,ly,my, el,ss, a)    ! [rbetat]
!=============================================================================
! Vector version of rbeta2t filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                                       intent(in   ):: nv, &
                                                               hx,Lx,mx, &
                                                               hy,ly,my
real(dp),dimension(   2,2,Lx:Mx,Ly:My),        intent(in   ):: el
real(dp),dimension(       Lx:Mx,Ly:My),        intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                            :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy):: b
real(dp),dimension(nv)                        :: ta
real(dp),dimension(2,2)                       :: tel
real(dp)                                      :: s,rr,rrx,rrc, &
                                                 frow,exx,eyy,eyx,x,y,xc
integer                                       :: ix,jx,gx
integer                                       :: iy,jy,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   ta=a(:,ix,iy); s=ss(ix,iy)
   tel=el(:,:,ix,iy)*this%rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy; y=gy; xc=-y*eyx
      rrx=(y*eyy)**2; x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**this%p
         b(:,jx,jy)=b(:,jx,jy)+frow*ta
      enddo! gx
   enddo! gy
enddo; enddo ! ix, iy
a=b
end subroutine vrbeta2t

!=============================================================================
module subroutine vrbeta3T(this,nv,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a) ! [rbetat]
!=============================================================================
! Vector version of rbeta3t filtering nv fields at once.
!=============================================================================
class(mg_parameter_type)::this
integer,                                    intent(in   ):: nv,      &
                                                            hx,Lx,mx,&
                                                            hy,ly,my,&
                                                            hz,lz,mz
real(dp),dimension(   3,3,Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: el
real(dp),dimension(       Lx:Mx,Ly:My,Lz:Mz),intent(in   ):: ss
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
                                lz-hz:mz+hz),intent(inout):: a
!-----------------------------------------------------------------------------
real(dp),parameter                         :: eps=1.e-12
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
                               lz-hz:mz+hz):: b
real(dp),dimension(nv)                     :: ta
real(dp),dimension(3,3)                    :: tel
real(dp):: s,rr,rrx,rry,rrc,frow,&
           exx,eyy,ezz,eyx,ezx,ezy,x,y,z,xc,yc
integer :: ix,jx,gx
integer :: iy,jy,gy
integer :: iz,jz,gz
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   ta=a(:,ix,iy,iz); s=ss(ix,iy,iz)
   tel=el(:,:,ix,iy,iz)*this%rmom2_3
   exx=tel(1,1); eyy=tel(2,2); ezz=tel(3,3)
   eyx=tel(2,1); ezx=tel(3,1); ezy=tel(3,2)
   z=u1/ezz
   do gz=ceiling(-z+eps),floor( z-eps)
      jz=iz+gz; z=gz; yc=-z*ezy
      rry=(z*ezz)**2; y =sqrt(u1-rry)
      do gy=ceiling((yc-y)/eyy+eps),floor((yc+y)/eyy-eps)
         jy=iy+gy; y=gy;        xc=-y*eyx-z*ezx
         rrx=rry+(y*eyy-yc)**2; x =sqrt(u1-rrx)
         do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
            jx=ix+gx; x=gx
            rr=rrx+(x*exx-xc)**2; rrc=u1-rr
            frow=s*rrc**this%p
            b(:,jx,jy,jz)=b(:,jx,jy,jz)+frow*ta
         enddo! gx
      enddo! gy
   enddo! gz
enddo; enddo; enddo! ix, iy, iz
a=b
end subroutine vrbeta3t

end submodule jp_pbfil
