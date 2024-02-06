!                                            *********************************
!                                            *      module pbfil
!                                            *
!                                            *      R. J. Purser
!                                            *
!                                            *      NOAA/NCEP/EMC
!                                            *
!                                            *      March 2019
!                                            *
!                                            *
!                                            *
!                                            *********************************
!
! Codes for the beta filters.
! The filters invoke the aspect tensor information encoded by the 
! Cholesky lower-triangular factors, el, of the INVERSE aspect tensors.
! The routines, "cholaspect", convert (in place) the field of given
! aspect tensors A to the equivalent cholesky factors of A^(-1).
! The routines, "getlinesum" precompute the normalization coefficients
! for
! each line (row) of the implied matrix form of the beta filter so that
! the
! normalized line sum associated with each point of application becomes
! unity. This makes the application of each filter significantly faster
! than having to work out the normalization on the fly.
!
! Be sure to have run cholaspect, and then getlinesum, prior to applying
! the
! beta filters themselves.
!
! Direct dependencies:
! Libraries: jp_pmat
! Modules:   jp_pkind, jp_pietc, jp_pmat
!            mg_parameter
!
!=============================================================================
module jp_pbfil
!=============================================================================
use mpi
use kinds, only: dp=>r_kind
!!!use jp_pkind, only: dp
use jp_pietc, only: u1
use mg_parameter, only: p, rmom2_1,rmom2_2,rmom2_3,rmom2_4
implicit none

private 

public cholaspect
interface cholaspect
   module procedure cholaspect1,cholaspect2,cholaspect3,cholaspect4
end interface

public getlinesum
interface getlinesum
   module procedure getlinesum1,getlinesum2,getlinesum3
end interface

public rbeta
interface rbeta
   module procedure rbeta1, rbeta2, rbeta3, rbeta4,    &
                   vrbeta1,vrbeta2,vrbeta3,vrbeta4                 
end interface

public rbetaT
interface rbetaT
   module procedure rbeta1t, rbeta2t, rbeta3t, rbeta4t, &
                   vrbeta1t,vrbeta2t,vrbeta3t,vrbeta4t
end interface


contains

!=============================================================================
subroutine cholaspect1(lx,mx, el)                              !  [cholaspect]
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
do ix=lx,mx; el(1,1,ix)=u1/sqrt(el(1,1,ix)); enddo
end subroutine cholaspect1
!=============================================================================
subroutine cholaspect2(lx,mx, ly,my, el)                        !  [cholaspect]
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
subroutine cholaspect3(lx,mx, ly,my, lz,mz, el)                 !  [cholaspect]
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
subroutine cholaspect4(lx,mx, ly,my, lz,mz, lw,mw,el)           !  [cholaspect]
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
subroutine getlinesum1(hx,lx,mx, el, ss)                       !  [getlinesum]
!=============================================================================
! Get inverse of the line-sum of the matrix representing the
! unnormalized
! beta function with aspect tensor pasp=(el*el^T)^(-1), and invert the
! result 
! so it can be used subsequently in the normalized version of this
! filter.
!=============================================================================
integer,                  intent(in   ):: hx,Lx,mx
real(dp),dimension(1,1,Lx:Mx),intent(in   ):: el
real(dp),dimension(lx:mx),intent(  out):: ss
!-----------------------------------------------------------------------------
real(dp),parameter:: eps=1.e-12
real(dp)          :: s,rr,rrc,exx,x
integer           :: ix,gxl,gxm,gx
!=============================================================================
do ix=Lx,Mx
   s=0
   exx=el(1,1,ix)*rmom2_1
   x=u1/exx
   gxl=ceiling(-x+eps); gxm=floor( x-eps)
   if(gxl<-hx.or.gxm>hx)&
        stop 'In getlinesum1; filter reach fx becomes too large for hx'
   do gx=gxl,gxm
      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      s=s+rrc**p
   enddo
   ss(ix)=u1/s
enddo
end subroutine getlinesum1
!=============================================================================
subroutine getlinesum2(hx,lx,mx, hy,ly,my, el, ss)             !  [getlinesum]
!=============================================================================
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
   tel=el(:,:,ix,iy)*rmom2_2 ! This el, rescaled
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
         s=s+rrc**p
      enddo! gx
   enddo! gy
   ss(ix,iy)=u1/s
enddo;  enddo! ix, iy
end subroutine getlinesum2
!=============================================================================
subroutine getlinesum3(hx,lx,mx, hy,ly,my, hz,lz,mz, el, ss)!     [getlinesum]
!=============================================================================
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
   tel=el(:,:,ix,iy,iz)*rmom2_3
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
            s=s+rrc**p
         enddo! gx
      enddo! gy
   enddo! gz
   ss(ix,iy,iz)=u1/s
enddo; enddo; enddo! ix, iy, iz
end subroutine getlinesum3
!=============================================================================
subroutine getlinesum4(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el, ss)                                                   !  [getlinesum]
!=============================================================================
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
   tel=el(:,:,ix,iy,iz,iw)*rmom2_4
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
               s=s+rrc**p
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
   ss(ix,iy,iz,iw)=u1/s
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
end subroutine getlinesum4

!=============================================================================
subroutine rbeta1(hx,lx,mx, el,ss, a)                               !  [rbeta]
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
   exx=el(ix)*rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**p
      tb=tb+frow*a(jx)
   enddo
   b(ix)=tb
enddo
a=b
end subroutine rbeta1
!=============================================================================
subroutine rbeta2(hx,lx,mx, hy,ly,my, el,ss, a)                     !  [rbeta]
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
   tel=el(:,:,ix,iy)*rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy;       y=gy; xc=-y*eyx
      rrx=(y*eyy)**2;       x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**p
         tb=tb+frow*a(jx,jy)
      enddo! gx
   enddo! gy
   b(ix,iy)=tb
enddo; enddo! ix, iy
a=b
end subroutine rbeta2
!=============================================================================
subroutine rbeta3(hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)            !  [rbeta]
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
   tel=el(:,:,ix,iy,iz)*rmom2_3
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
            frow=s*rrc**p
            tb=tb+frow*a(jx,jy,jz)
         enddo! gx
      enddo! gy
   enddo! gz
   b(ix,iy,iz)=tb
enddo;   enddo;    enddo! ix, iy, iz
a=b
end subroutine rbeta3
!=============================================================================
subroutine rbeta4(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,ss,a)!    [rbeta]
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
   tel=el(:,:,ix,iy,iz,iw)*rmom2_4
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
               frow=s*rrc**p
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
subroutine vrbeta4(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss,a)                                                       !  [rbeta]
!=============================================================================
! Vector version of rbeta4 filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy,iz,iw)*rmom2_4
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
               frow=s*rrc**p
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
subroutine rbeta1T(hx,lx,mx, el,ss, a)                             !  [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 1D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= mx+hx.
!=============================================================================
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
   exx=el(1,1,ix)*rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**p
      b(jx)=b(jx)+frow*ta
   enddo
enddo
a=b
end subroutine rbeta1t
!=============================================================================
subroutine rbeta2T(hx,lx,mx, hy,ly,my, el,ss, a)                   !  [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 2D.
! It conserved "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx, Ly <= iy <= My.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= mx+hx, Ly-hy <= Jy <= my+hy
!=============================================================================
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
   tel=el(:,:,ix,iy)*rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy; y=gy; xc=-y*eyx
      rrx=(y*eyy)**2; x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**p
         b(jx,jy)=b(jx,jy)+frow*ta
      enddo! gx 
   enddo! gy
enddo;  enddo! ix, iy
a=b
end subroutine rbeta2t
!=============================================================================
subroutine rbeta3T(hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)         !  [rbetat]
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 3D.
! It conserves "masses" initially distributed only at the closure of 
! the central domain, 
! Lx <= ix <= Mx, Ly <= iy <= My, Lz <= iz <= Mz.
! The output field of the redistributed masses occupies the
! the extended domain, 
! Lx-hx <= jx <= Mx+hx, Ly-hy <= Jy <= My+hy, Lz-hz <= Jz <= Mz+hz.
!=============================================================================
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
   tel=el(:,:,ix,iy,iz)*rmom2_3
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
            frow=s*rrc**p
            b(jx,jy,jz)=b(jx,jy,jz)+frow*ta
         enddo! gx
      enddo! gy
   enddo ! gz
enddo;  enddo;  enddo ! ix, iy, iz
a=b
end subroutine rbeta3t
!=============================================================================
subroutine rbeta4T(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     el,ss, a)                                                   !  [rbetat]
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
   tel=el(:,:,ix,iy,iz,iw)*rmom2_4
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
               frow=s*rrc**p
               b(jx,jy,jz,jw)=b(jx,jy,jz,jw)+frow*ta
            enddo! gx
         enddo! gy
      enddo! gz
   enddo! gw
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
a=b
end subroutine rbeta4t


!=============================================================================
subroutine vrbeta4t(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, &
                                              hw,lw,mw, el,ss, a)!    [rbetat]
!=============================================================================
! Vector version of rbeta4t filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy,iz,iw)*rmom2_4
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
               frow=s*rrc**p
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
subroutine vrbeta1(nv,hx,lx,mx, el,ss, a)                           !  [rbeta]
!=============================================================================
! Vector version of rbeta1 filtering nv fields at once.
!=============================================================================
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
   exx=el(1,1,ix)*rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**p
      tb=tb+frow*a(:,jx)
   enddo
   b(:,ix)=tb
enddo
a=b
end subroutine vrbeta1

!=============================================================================
subroutine vrbeta2(nv,hx,lx,mx, hy,ly,my, el,ss, a)                 !  [rbeta]
!=============================================================================
! Vector version of rbeta2 filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy)*rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy;       y=gy; xc=-y*eyx
      rrx=(y*eyy)**2;       x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**p
         tb=tb+frow*a(:,jx,jy)
      enddo! gx
   enddo! gy
   b(:,ix,iy)=tb
enddo;   enddo! ix, iy
a=b
end subroutine vrbeta2

subroutine vrbeta3(nv, hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss,a)       !  [rbeta]
!=============================================================================
! Vector version of rbeta3 filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy,iz)*rmom2_3
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
            frow=s*rrc**p
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
subroutine vrbeta1T(nv, hx,lx,mx, el,ss, a)                        !  [rbetat]
!=============================================================================
! Vector version of rbeta1t filtering nv fields at once.
!=============================================================================
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
   exx=el(1,1,ix)*rmom2_1
   x=u1/exx
   do gx=ceiling(-x+eps),floor( x-eps)
      jx=ix+gx;      x=gx
      rr=(x*exx)**2; rrc=u1-rr
      frow=s*rrc**p
      b(:,jx)=b(:,jx)+frow*ta
   enddo
enddo
a=b
end subroutine vrbeta1t
!=============================================================================
subroutine vrbeta2T(nv,hx,lx,mx, hy,ly,my, el,ss, a)               !  [rbetat]
!=============================================================================
! Vector version of rbeta2t filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy)*rmom2_2 ! This el, rescaled
   exx=tel(1,1); eyy=tel(2,2)
   eyx=tel(2,1)
   y=u1/eyy
   do gy=ceiling(-y+eps),floor( y-eps)
      jy=iy+gy; y=gy; xc=-y*eyx
      rrx=(y*eyy)**2; x =sqrt(u1-rrx)
      do gx=ceiling((xc-x)/exx+eps),floor((xc+x)/exx-eps)
         jx=ix+gx; x=gx
         rr=rrx+(x*exx-xc)**2; rrc=u1-rr
         frow=s*rrc**p
         b(:,jx,jy)=b(:,jx,jy)+frow*ta
      enddo! gx
   enddo! gy
enddo; enddo ! ix, iy
a=b
end subroutine vrbeta2t

!=============================================================================
subroutine vrbeta3T(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, el,ss, a)     !  [rbetat]
!=============================================================================
! Vector version of rbeta3t filtering nv fields at once.
!=============================================================================
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
   tel=el(:,:,ix,iy,iz)*rmom2_3
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
            frow=s*rrc**p
            b(:,jx,jy,jz)=b(:,jx,jy,jz)+frow*ta
         enddo! gx
      enddo! gy
   enddo! gz
enddo; enddo; enddo! ix, iy, iz
a=b
end subroutine vrbeta3t

end module jp_pbfil

