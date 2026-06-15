!#
!                                            *********************************
!                                            *      module wbfil             *
!                                            *      R. J. Purser             *
!                                            *      NOAA/NCEP/EMC            *
!                                            *      May 2026                 *
!                                            *    jim.purser@noaa.gov        *
!                                            *********************************
!
! Codes for the beta line and radial filters.
! The filters invoke the aspect tensor information encoded by the 
! Cholesky lower-triangular factors, el, of the INVERSE aspect tensors.
! The routines, "rcalib", convert the field, as, of given
! aspect tensors A to the equivalent cholesky factors of A^(-1),
! and get the normalization coefficients for
! each line (row) of the implied matrix form of the beta filter so that the
! normalized line sum associated with each point of application becomes
! unity. This makes the application of each filter significantly faster
! than having to work out the normalization on the fly.
! This is version adapted from ybfil.f90, but with a more general specification 
! of boundary normal aspect tensor information for use in the cosmetic boundary
! treatment in limited domains. Also, this version should include the vertical
! multigrid filter option.
!
! COMPILE AFTER: {pmat}
!
!=============================================================================
module wbfil
!=============================================================================
use pkind, only: spi,dp
use pietc, only: u0,u1
implicit none
private
public:: &
     t22_to_3,t2_to_3,t3_to_22,t33_to_6,t3_to_6,t6_to_33,&
     t44_to_10,t4_to_10,t10_to_44,p,rsg,i2pair,i3pair,i4pair, &
     finp,inip,&
     rcalib,rbeta,rbetat
integer(spi),parameter       :: nbcof=9,np=6
real(dp),parameter           :: eps=1.e-12_dp,emin=.3e0_dp
real(dp),dimension(0:nbcof,0:nbcof):: bcofs
real(dp)                     :: op,rsg,rpp3o2,rpp4o2,rpp5o2,rpp6o2,nm2x
integer(spi)                 :: p,p2p3,p2p4
integer(spi),dimension(2,0:2):: i2pair
integer(spi),dimension(2,6)  :: i3pair
integer(spi),dimension(2,10) :: i4pair
data p/0/ ! default (the filter will fail with p=0, alerting the user to call inip)
data i2pair/1,1, 2,2, 2,1/
data i3pair/1,1, 2,2, 3,3, 3,2, 3,1, 2,1/
data i4pair/1,1, 2,2, 3,3, 4,4, 2,1, 3,1, 4,1, 3,2, 4,2, 4,3/
interface t22_to_3;    module procedure i22_to_3, r22_to_3;   end interface
interface t2_to_3;     module procedure i2_to_3,  r2_to_3;    end interface
interface t3_to_22;    module procedure i3_to_22, r3_to_22;   end interface
interface t33_to_6;    module procedure i33_to_6, r33_to_6;   end interface
interface t3_to_6;     module procedure i3_to_6,  r3_to_6;    end interface
interface t6_to_33;    module procedure i6_to_33, r6_to_33;   end interface
interface t44_to_10;   module procedure i44_to_10,r44_to_10;  end interface
interface t4_to_10;    module procedure i4_to_10, r4_to_10;   end interface
interface t10_to_44;   module procedure i10_to_44,r10_to_44;  end interface
!---
interface finp;        module procedure finp;                 end interface
interface inip;        module procedure inip;                 end interface
interface hofnm2;      module procedure hofnm2;               end interface
interface nm2ofh;      module procedure nm2ofh;               end interface
interface bfmoms;      module procedure bfmoms;               end interface
interface rcalib;      module procedure rcalib1,rcalib2;      end interface
interface rbetat
   module procedure rbeta1t,vrbeta1t, rbeta2t,vrbeta2t
end interface
interface rbeta
   module procedure rbeta1,vrbeta1,   rbeta2,vrbeta2
end interface
contains

!==============================================================================
subroutine i22_to_3(i22,i3)!                                         [t22_to_3]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(2,2),intent(in ):: i22
integer(spi),dimension(0:2),intent(out):: i3
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2; i3(L)=i22(i2pair(1,L),i2pair(2,L)); enddo
end subroutine i22_to_3
!==============================================================================
subroutine r22_to_3(r22,r3)!                                         [t22_to_3]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(2,2),intent(in ):: r22
real(dp),dimension(0:2),intent(out):: r3
!-----------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2; r3(L)=r22(i2pair(1,L),i2pair(2,L)); enddo
end subroutine r22_to_3

!==============================================================================
subroutine i2_to_3(i2,i3)!                                            [t2_to_3]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(2),intent(in ):: i2
integer(spi),dimension(3),intent(out):: i3
!------------------------------------------------------------------------------
call t22_to_3(outer_product(i2,i2),i3)
end subroutine i2_to_3
!==============================================================================
subroutine r2_to_3(r2,r3)!                                            [t2_to_3]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(2),intent(in ):: r2
real(dp),dimension(3),intent(out):: r3
!------------------------------------------------------------------------------
call t22_to_3(outer_product(r2,r2),r3)
end subroutine r2_to_3

!==============================================================================
subroutine i3_to_22(i3,i22)!                                         [t3_to_22]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(0:2),intent(in ):: i3
integer(spi),dimension(2,2),intent(out):: i22
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2
   i22(i2pair(1,L),i2pair(2,L))=i3(L)
   i22(i2pair(2,L),i2pair(1,L))=i3(L)
enddo
end subroutine i3_to_22
!==============================================================================
subroutine r3_to_22(r3,r22)!                                         [t3_to_22]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(0:2),intent(in ):: r3
real(dp),dimension(2,2),intent(out):: r22
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2
   r22(i2pair(1,L),i2pair(2,L))=r3(L)
   r22(i2pair(2,L),i2pair(1,L))=r3(L)
enddo
end subroutine r3_to_22

!==============================================================================
subroutine i33_to_6(i33,i6)!                                         [t33_to_6]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(3,3),intent(in ):: i33
integer(spi),dimension(6)  ,intent(out):: i6
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6; i6(L)=i33(i3pair(1,L),i3pair(2,L)); enddo
end subroutine i33_to_6
!==============================================================================
subroutine r33_to_6(r33,r6)!                                         [t33_to_6]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(3,3),intent(in ):: r33
real(dp),dimension(6)  ,intent(out):: r6
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6; r6(L)=r33(i3pair(1,L),i3pair(2,L)); enddo
end subroutine r33_to_6

!==============================================================================
subroutine i3_to_6(i3,i6)!                                            [t3_to_6]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(3),intent(in ):: i3
integer(spi),dimension(6),intent(out):: i6
!------------------------------------------------------------------------------
call t33_to_6(outer_product(i3,i3),i6)
end subroutine i3_to_6
!==============================================================================
subroutine r3_to_6(r3,r6)!                                            [t3_to_6]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(3),intent(in ):: r3
real(dp),dimension(6),intent(out):: r6
!------------------------------------------------------------------------------
call t33_to_6(outer_product(r3,r3),r6)
end subroutine r3_to_6

!==============================================================================
subroutine i6_to_33(i6,i33)!                                         [t6_to_33]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(6),  intent(in ):: i6
integer(spi),dimension(3,3),intent(out):: i33
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6
   i33(i3pair(1,L),i3pair(2,L))=i6(L)
   i33(i3pair(2,L),i3pair(1,L))=i6(L)
enddo
end subroutine i6_to_33
!==============================================================================
subroutine r6_to_33(r6,r33)!                                         [t6_to_33]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(6),  intent(in ):: r6
real(dp),dimension(3,3),intent(out):: r33
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6
   r33(i3pair(1,L),i3pair(2,L))=r6(L)
   r33(i3pair(2,L),i3pair(1,L))=r6(L)
enddo
end subroutine r6_to_33

!==============================================================================
subroutine i44_to_10(i44,i10)!                                      [t44_to_10]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(4,4),intent(in ):: i44
integer(spi),dimension(10) ,intent(out):: i10
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10; i10(L)=i44(i4pair(1,L),i4pair(2,L)); enddo
end subroutine i44_to_10
!==============================================================================
subroutine r44_to_10(r44,r10)!                                      [t44_to_10]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(4,4),intent(in ):: r44
real(dp),dimension(10) ,intent(out):: r10
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10; r10(L)=r44(i4pair(1,L),i4pair(2,L)); enddo
end subroutine r44_to_10

!==============================================================================
subroutine i4_to_10(i4,i10)!                                         [t4_to_10]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(4), intent(in ):: i4
integer(spi),dimension(10),intent(out):: i10
!------------------------------------------------------------------------------
call t44_to_10(outer_product(i4,i4),i10)
end subroutine i4_to_10
!==============================================================================
subroutine r4_to_10(r4,r10)!                                         [t4_to_10]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(4), intent(in ):: r4
real(dp),dimension(10),intent(out):: r10
!------------------------------------------------------------------------------
call t44_to_10(outer_product(r4,r4),r10)
end subroutine r4_to_10

!==============================================================================
subroutine i10_to_44(i10,i44)!                                      [t10_to_44]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(10), intent(in ):: i10
integer(spi),dimension(4,4),intent(out):: i44
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10
   i44(i4pair(1,L),i4pair(2,L))=i10(L)
   i44(i4pair(2,L),i4pair(1,L))=i10(L)
enddo
end subroutine i10_to_44
!==============================================================================
subroutine r10_to_44(r10,r44)!                                      [t10_to_44]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(10), intent(in ):: r10
real(dp),dimension(4,4),intent(out):: r44
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10
   r44(i4pair(1,L),i4pair(2,L))=r10(L)
   r44(i4pair(2,L),i4pair(1,L))=r10(L)
enddo
end subroutine r10_to_44

!--------------------------------------------------------

!======================================================================= [finp]
subroutine finp
!==============================================================================
implicit none
p=0 ! revert to default value
end subroutine finp
!======================================================================= [inip]
subroutine inip(p_prescribe,ff)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u2,u3,u4,u5,o2,o3
implicit none
integer(spi),intent(in ):: p_prescribe
logical,     intent(out):: ff
!------------------------------------------------------------------------------
real(dp),parameter                     :: u16=16,u4o3=u4*o3,u3o2=u3*o2,u5o2=u5*o2
real(dp),dimension(np*4+1)             :: ffac
real(dp),dimension(np)                 :: fac
integer(spi),dimension(0:nbcof,0:nbcof):: inums! numerators for B coefficients
integer(spi),dimension(0:nbcof)        :: idens! denominators for B coefficients
integer(spi)                           :: i
data inums/2,9*0, 1,2,8*0, -1,10,6,7*0, 1,-7,21,6,6*0, -3,20,-42,60,10,5*0,&
     5,-33,66,-66,55,6,4*0, -691,4550,-9009,8580,-5005,2730,210,3*0,       &
     105,-691,1365,-1287,715,-273,105,6,0,0,                               &
     -3617,23800,-46988,44200,-24310,8840,-2380,680,30,0,                  &
     219335,-1443183,2848860,-2678316,1469650,-529074,135660,-27132,5985,210/
data idens/1,3,15,21,45,33,1365,45,255,1995/
!==============================================================================
ff=(p_prescribe<1 .or. p_prescribe>np)
if(ff)then
   print'(" In inip; prescribed exponent p out of bounds")'
   return
endif
p=p_prescribe
op=u1/p
nm2x=u2/(u4o3**p+u2)! The normalized 2nd moment for a line-filter of half-span of 2.
p2p3=p*2+3
p2p4=p*2+4
rpp3o2=sqrt(p+u3o2)
rpp4o2=sqrt(p+u2)
rpp5o2=sqrt(p+u5o2)
rpp6o2=sqrt(p+u3)
ffac(1)=1
do i=3,p*4+1,2
   ffac(i)=ffac(i-2)*i
enddo
fac(1)=1
do i=2,p
   fac(i)=fac(i-1)*i
enddo
rsg=ffac(p*4+1)*u2**(p-1)/(ffac(p*2-1)*fac(p)*u16**p)
bcofs=inums; do i=0,nbcof; bcofs(:,i)=bcofs(:,i)/idens(i); enddo

end subroutine inip

!======================================================================[hofnm2]
subroutine hofnm2(nm2t,h)
!===============================================================================
! For nm2t<nm2x, the half-span h is less than 2 and h can then be computed directly.
! Otherwise we use Newton iterations to seek the half-span, h, such that the single
! application of the beta line filter with the integer exponent, p, previously set
! in the initialization routine inip produces the target normalized 2nd moment, nm2t.
!===============================================================================
use pkind, only: dp
use pietc, only: u1
implicit none
real(dp),intent(in ):: nm2t! target normalized 2nd-moment
real(dp),intent(out):: h   ! half-span of beta line filter
!-------------------------------------------------------------------------------
integer,parameter :: nit=40! Maximum limit for the number of Newton iterations
real(dp)          :: nm2,dnm2,r
integer           :: it
!===============================================================================
if(nm2t<nm2x)then ! Qualifying for the special cases where the half span h <= 2:
   h=u1/sqrt(u1-(nm2t/(2*(u1-nm2t)))**op)! <- analytic formula in special cases
   return
endif
! For the cases where the half-span will be found to be h >= 2: 
! p2p3 is the asymptotic limit of h**2/nm2 as h and nm2 go to infinity:
h=sqrt(u1+p2p3*nm2t)! <- First guess of half-span initializing Newton iterations
do it=1,nit
   call nm2ofh(h,nm2,dnm2)! Get normalized 2nd moment,nm2, for half-span, h.
   r=(nm2-nm2t)/dnm2! <- residual error of nm2 divided by d(nm2)/dh
   h=h-r
   if(abs(r)<eps)return ! <- return when the last correction was small enough
enddo
end subroutine hofnm2

!======================================================================[nm2ofh]
subroutine nm2ofh(h,nm2,dnm2)
!=============================================================================
! Compute the exact normalize 2nd moment, nm2, and its derivative dnm2,
! of the beta line filter of half-span h on a unit grid when the filter
! exponent parameter is the p initialized in subr. inip.
!=============================================================================
use pkind, only: dp
use pietc, only: u1
implicit none
real(dp),intent(in ):: h
real(dp),intent(out):: nm2,dnm2
!-----------------------------------------------------------------------------
real(dp):: mom0,mom2,dmom0,dmom2,omom0
!=============================================================================
call bfmoms(h,mom0,mom2,dmom0,dmom2)
omom0=u1/mom0
nm2=mom2*omom0
dnm2=(dmom2*mom0-dmom0*mom2)*omom0**2
end subroutine nm2ofh

!======================================================================[bfmoms]
subroutine bfmoms(h,mom0,mom2,dmom0,dmom2)
!==============================================================================
! For a beta line filter of exponent p and half-span h, compute the
! exact 0th and 2nd moments, mom0 and mom2, together with their
! derivatives wrt h, dmom0 and dmom2.
! The method recognizes the beta filter profile to be an even-polynomial
! so that a residual-free Euler-Maclaurin expansion can be used to
! evaluate the discrete summation in terms of the corresponding integral
! and the finitely many end-correction terms. The sum of each even-power
! of a symmetric range of integers, [-n:n] is that even power of n plus
! an odd polynomial in n. The coefficients of that odd polynomial are
! tabulated, for each of the original even powers, in array bcofs. 
!==============================================================================
use pkind, only: dp
use pietc, only: u1,o2
implicit none
real(dp),intent(in ):: h
real(dp),intent(out):: mom0,mom2,dmom0,dmom2
!------------------------------------------------------------------------------
real(dp),dimension(0:np)     :: pchoose
real(dp)                     :: ho2,hh,c,dc,q0,q2,enn
integer                      :: i,j,jm,jp,k,n
!==============================================================================
! Set up binomial coefficients p-choose-j:
pchoose(0)=u1
do j=1,p; jm=j-1; pchoose(j)=(pchoose(jm)*(p-jm))/j; enddo
ho2=h*o2; hh=h*h; n=h; enn=n*n
mom0=0; dmom0=0; mom2=0; dmom2=0
do j=0,p; jp=j+1
   c=pchoose(j)/(-hh)**j; dc=-j*c/ho2
   q0=enn**j;             q2=q0*enn
   do k=0,j;  q0=q0+n*bcofs(k,j )*enn**k; enddo! Untruncated Euler-Maclaurin formula needs no
   do k=0,jp; q2=q2+n*bcofs(k,jp)*enn**k; enddo! residual and becomes finite and exact for even powers.
   mom0  =mom0 +c*q0;    mom2= mom2+ c*q2
   dmom0=dmom0+dc*q0;   dmom2=dmom2+dc*q2
enddo
end subroutine bfmoms

! Calibration routines to initialize filter scale constants at each treated boundary,
!  the amplitude normalization coefficients, el(0,..), the Cholesky coefficients for
! each given overall aspect tensor, and the halo widths needed by the filters in
! each particular subdomain.
!=========================================================================[rcalib]
subroutine rcalib1(hx,Lx,mx,as,el,hxm)!                            
!=================================================================================
! For each given aspect tensor, as, find the exact half-span such that the
! single application of the beta filter of this half-span produces the normalized
! second moment, as/2. This implies that the adjoint, followed by the direct,
! beta filter will result in the desired aspect tensor, as, overall.
! The computed half-span must be strictly less than the integer, hx+1, where
! hx is the prescribed halo width in the direction of filtering.
!
! The reciprocal of each computed half span is stored in EL(1,:), while the
! computed normalizing coefficient, taken from the square-root 
! of the inferred adjoint+direct unnormalized filter amplitude. is stored
! in EL(0,:).
!
! For true non-oblique line filters, it is never necessary to apply them to
! points in the halos, so the EL coefficents are only evluated in the interior
! of the subdomain.
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u1,u2,o2
use pmat,  only: inv, L1Lm
implicit none
integer(spi),                 intent(in ):: hx,Lx,mx
real(dp),dimension(    Lx:mx),intent(in ):: as
real(dp),dimension(0:1,Lx:mx),intent(out):: eL
integer(spi),dimension(Lx:mx),intent(out):: hxm
!------------------------------------------------------------------------------
real(dp),dimension(-hx:hx):: fs
real(dp)                  :: exx,f,r,rc,rrc,s
integer(spi)              :: ix,gx,gxm,gxn
!==============================================================================
if(p==0)stop 'In rcalib1; Error: Initialization routine, inip, has not been called to initialize p'
do ix=Lx,mx
   call hofnm2(as(ix)*o2,s)
   exx=u1/s
   el(1,ix)=exx
   gxm=floor(u1/exx); hxm(ix)=gxm
   fs(-gxm:gxm)=0
   fs(0)=u1
   do gx=-gxm,-1
      rrc=u1-(gx*exx)**2
      f=rrc**p; fs(-gx)=f; fs(gx)=f
   enddo
   eL(0,ix)=u1/sqrt(sum(fs(-gxm:gxm)**2))
enddo
end subroutine rcalib1

!===================================================================== [rcalib]
subroutine rcalib2(hx,Lx,mx, hy,Ly,my, as, eL, hxm,hym)
!==============================================================================
! Convert the given field, as, of aspect tensors into the equivalent field el(1:,..)
! of Cholesky lower-triangular factors of the inverses of the aspect tensors
! in 2D, and put the exact normalizing coefficient in eL(0,..)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u2,o2
use pmat, only: inv, L1Lm
implicit none
integer(spi),                       intent(in ):: hx,Lx,mx, hy,Ly,my
real(dp),dimension(3,  Lx:mx,Ly:my),intent(in ):: as
real(dp),dimension(0:3,Lx:mx,Ly:my),intent(out):: eL
integer(spi),dimension(Lx:mx,Ly:my),intent(out):: hxm,hym
!------------------------------------------------------------------------------
real(dp),parameter:: epss=1.e-30
real(dp),dimension(-hx:hx,-hy:hy):: fs
real(dp),dimension(2,2):: tas,tel22
real(dp),dimension(3)  :: tel
real(dp)               :: aa0,cx,exx,eyy,eyx,f,r,rrc,rrxc
integer(spi)           :: gx,gxL,gxm,gxn, gy,gyL,gym,gyn, &
     gxmm,&
     ix,ixp,ixm, iy,iyp,iym, &
     kx,nx,     ky,ny
!==============================================================================
if(p==0)stop 'In rcalib2; Error: Initialization routine, inip, has not been called to initialize p'
fs=0
do iy=Ly,my; do ix=Lx,mx
   call t3_to_22(as(:,ix,iy),tas)
   tas=p2p4*tas/2; call inv(tas); call L1Lm(tas,tel22)
   call t22_to_3(tel22,tel)
   el(1:3,ix,iy)=tel
   exx=tel(1); eyy=tel(2); eyx=tel(3)
   fs(0,0)=u1
   gym=floor(u1/eyy); gxmm=0
   lgy: do gy=-gym,0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      gxL=ceiling((cx-r)/exx); gxm=floor((cx+r)/exx); gxmm=max(-gxL,max(gxm,gxmm))
      do gx=gxL,gxm; ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         f=rrc**p
         fs(gx,gy)=f; fs(-gx,-gy)=f
      enddo! gx
   enddo lgy
   gxm=gxmm
   hxm(ix,iy)=gxm
   hym(ix,iy)=gym

   eL(0,ix,iy)=u1/sqrt(sum(fs(-gxm:gxm,-gym:gym)**2))
   fs(-gxm:gxm,-gym:gym)=0
enddo; enddo!  ix iy
end subroutine rcalib2

!# Adjoint beta filters:
!======================================================================[rbetat]
subroutine rbeta1T(hx,Lx,mx, el, a,b)
!=============================================================================
! Perform an ADJOINT beta-function filter in 1D.
!=============================================================================
integer(spi),                    intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:mx   ),intent(in   ):: el
real(dp),dimension(    Lx:mx   ),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp)    :: tafrow,tas
real(dp)    :: exx,rrc
integer(spi):: ix,ixp,ixm,gx
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tas=a(ix)*el(0,ix)
   b(ix)=b(ix)+tas
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tafrow=tas*rrc**p
      b(ixp)=b(ixp)+tafrow
      b(ixm)=b(ixm)+tafrow
   enddo
enddo
end subroutine rbeta1t
!======================================================================[rbetat]
subroutine vrbeta1T(nv, hx,lx,mx, el, a,b)
!=============================================================================
! Vector version of rbeta1t filtering nv fields at once.
!=============================================================================
integer(spi),                      intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(0:1,  Lx:mx),   intent(in   ):: el
real(dp),dimension(nv,   Lx:mx),   intent(in   ):: a
real(dp),dimension(nv,Lx-hx:mx+hx),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(nv):: tafrow,tas
real(dp)              :: exx,rrc
integer(spi)          :: ix,ixp,ixm,gx
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tas=a(:,ix)*el(0,ix)
   b(:,ix)=b(:,ix)+tas
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tafrow=tas*rrc**p
      b(:,ixp)=b(:,ixp)+tafrow
      b(:,ixm)=b(:,ixm)+tafrow
   enddo
enddo
end subroutine vrbeta1t
!======================================================================[rbetat]
subroutine rbeta2T(hx,lx,mx, hy,ly,my, el, a,b)
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 2D.
!=============================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my
real(dp),dimension(0:3,Lx:mx,       Ly:my   ),intent(in   ):: el
real(dp),dimension(    Lx:mx,       Ly:my   ),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3):: tel
real(dp)             :: tafrow,tas
real(dp)             :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)         :: gx,gy,ix,ixp,ixm,iy,iyp,iym
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tas=a(ix,iy)*el(0,ix,iy)
   b(ix,iy)=b(ix,iy)+tas
Lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit Lgy
         rrc=rrxc-(gx*exx-cx)**2
         tafrow=tas*rrc**p
         b(ixp,iyp)=b(ixp,iyp)+tafrow
         b(ixm,iym)=b(ixm,iym)+tafrow
      enddo! gx
   enddo Lgy
enddo;  enddo! ix, iy
end subroutine rbeta2t
!======================================================================[rbetat]
subroutine vrbeta2T(nv,hx,Lx,mx, hy,Ly,my, el, a,b)
!=============================================================================
! Vector version of rbeta2t filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,Ly,my
real(dp),dimension( 0:3,  Lx:mx,       Ly:my),   intent(in   ):: eL
real(dp),dimension(nv,    Lx:mx,       Ly:my),   intent(in   ):: a
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3) :: tel
real(dp),dimension(nv):: tafrow,tas
real(dp)              :: cx,rrc,rrxc,exx,eyy,eyx,r
integer(spi)          :: ix,ixp,ixm,gx,iy,iyp,iym,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tas=a(:,ix,iy)*el(0,ix,iy)
   b(:,ix,iy)=b(:,ix,iy)+tas
Lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit Lgy
         rrc=rrxc-(gx*exx-cx)**2
         tafrow=tas*rrc**p
         b(:,ixp,iyp)=b(:,ixp,iyp)+tafrow
         b(:,ixm,iym)=b(:,ixm,iym)+tafrow
      enddo! gx
   enddo Lgy
enddo; enddo ! ix, iy
end subroutine vrbeta2t

!# Direct beta filters:
!=========================================================================[rbeta]
subroutine rbeta1(hx,Lx,mx, el, a,b)
!===============================================================================
! Perform a direct beta-function filter in 1D.
!
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx.
! The output data occupy the central region
! Lx <= ix <= Mx.
!===============================================================================
use pkind, only: dp,spi
use pietc, only: u1
implicit none
integer(spi),                    intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:mx   ),intent(in   ):: el
real(dp),dimension(-hx+Lx:mx+hx),intent(in   ):: a
real(dp),dimension(    Lx:mx   ),intent(  out):: b
!-------------------------------------------------------------------------------
real(dp)     :: tb
real(dp)     :: exx,rrc
integer(spi) :: gx,ix,ixp,ixm
!===============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(ix)
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**p*(a(ixp)+a(ixm))
   enddo
   b(ix)=tb*el(0,ix)
enddo
end subroutine rbeta1
!===================================================================[rbeta]
subroutine vrbeta1(nv,hx,lx,mx, el, a,b)
!================================================================================
! Vector version of rbeta1 filtering nv fields at once.
!=============================================================================
integer(spi),                      intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(0:1,Lx:mx),     intent(in   ):: el
real(dp),dimension(nv,lx-hx:mx+hx),intent(in   ):: a
real(dp),dimension(nv, Lx:mx),     intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(nv):: tb
real(dp)              :: exx,rrc
integer(spi)          :: gx,ix,ixp,ixm
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(:,ix)
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**p*(a(:,ixp)+a(:,ixm))
   enddo
   b(:,ix)=tb*el(0,ix)
enddo
end subroutine vrbeta1
!========================================================================[rbeta]
subroutine rbeta2(hx,Lx,mx, hy,Ly,my, el, a,b)
!==============================================================================
! Perform a direct radial beta-function filter in 2D.
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u1
implicit none
integer(spi),                                 intent(in   ):: hx,Lx,mx,hy,Ly,my
real(dp),dimension(0:3,Lx:mx,       Ly:my   ),intent(in   ):: el
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(in   ):: a
real(dp),dimension(    Lx:mx,       Ly:my   ),intent(  out):: b
!------------------------------------------------------------------------------
real(dp),dimension(3):: tel
real(dp)             :: tb
real(dp)             :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)         :: gx,gy,ix,ixm,ixp,iy,iym,iyp
!==============================================================================
b=0
do iy=Ly,my; do ix=Lx,mx
   tel=el(1:3,ix,iy) 
   exx=tel(1); eyy=tel(2); eyx=tel(3)
   tb=a(ix,iy)
lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx);ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         tb=tb+rrc**p*(a(ixp,iyp)+a(ixm,iym))
      enddo! gx
   enddo lgy
   b(ix,iy)=tb*el(0,ix,iy)
enddo; enddo! ix, iy
end subroutine rbeta2
!====================================================================== [rbeta]
subroutine vrbeta2(nv,hx,lx,mx, hy,ly,my, el, a,b)
!=============================================================================
! Vector version of rbeta2 filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my
real(dp),dimension(0:3,Lx:Mx,Ly:My),           intent(in   ):: el
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(in   ):: a
real(dp),dimension(nv, lx:mx,ly:my),           intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3) :: tel
real(dp),dimension(nv):: tb
real(dp)              :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)          :: gx,gy,ix,ixp,ixm,iy,iyp,iym
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tb=a(:,ix,iy)
lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         tb=tb+rrc**p*(a(:,ixp,iyp)+a(:,ixm,iym))
      enddo! gx
   enddo lgy
   b(:,ix,iy)=tb*el(0,ix,iy)
enddo;   enddo! ix, iy
end subroutine vrbeta2

end module wbfil
!#
