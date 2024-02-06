!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_intstate)  mg_interpolate
!***********************************************************************
!                                                                      !
!    general mapping between 2d arrays using linerly squared           !
!    interpolations                                                    !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use kinds
use jp_pkind2, only: fpi

!use mpimod, only: mype
implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      module    subroutine def_offset_coef (this)                        
!***********************************************************************
implicit none
class(mg_intstate_type),target::this
                                                                      
real(r_kind):: r64,r32,r128,r2
!-----------------------------------------------------------------------
 r64 = 1.0d0/64.0d0
 r32 = 1.0d0/32.0d0
 r128= 1.0d0/128.0d0
! r2 = 1.0d0/2.0d0
 r2 = 1.0d0

! p_coef =(/-3.,51,29,-3/)
! q_coef =(/-3.,19.0d0,51.0d0,-3.0d0/)
! p_coef = p_coef*r64
! q_coef = q_coef*r64

 this%p_coef =(/-9.,111.,29.,-3./)
 this%q_coef =(/-3.,29.,111.,-9./)
this%p_coef = this%p_coef*r128 *r2
 this%q_coef = this%q_coef*r128 *r2

 this%a_coef =(/5.0d0,30.0d0,-3.0d0/)
 this%b_coef =(/-3.0d0,30.0d0,5.0d0/)
 this%a_coef=this%a_coef*r32 *r2
 this%b_coef=this%b_coef*r32 *r2
!-----------------------------------------------------------------------
                        endsubroutine def_offset_coef

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             module           subroutine lsqr_mg_coef (this)                           
!***********************************************************************
!                                                                      !
!   Prepare coeficients for mapping between:                           !
!   filter grid on analysis decomposition:  W(i0-ib:im+ib,j0-jb:jm+jb) !
!   and analysis grid:                      V(1:nm,1:mm)               !  
!                       - offset version -                             !
!                                                                      !
!              (  im < nm  and  jm < mm   )                            !
!                                                                      !
!***********************************************************************
implicit none
class(mg_intstate_type),target::this
real(r_kind), dimension(1:this%nm):: xa
real(r_kind), dimension(1-this%ib:this%im+this%ib):: xf
real(r_kind), dimension(1:this%mm):: ya
real(r_kind), dimension(1-this%jb:this%jm+this%jb):: yf
integer(i_kind):: i,j,n,m
real(r_kind) x1,x2,x3,x4,x
real(r_kind) x1x,x2x,x3x,x4x
real(r_kind) rx2x1,rx3x1,rx4x1,rx3x2,rx4x2,rx4x3
real(r_kind) y1,y2,y3,y4,y
real(r_kind) y1y,y2y,y3y,y4y
real(r_kind) ry2y1,ry3y1,ry4y1,ry3y2,ry4y2,ry4y3
real(r_kind) cfl1,cfl2,cfl3,cll
real(r_kind) cfr1,cfr2,cfr3,crr
!-----------------------------------------------------------------------
!
! Initialize
!
 
   do n=1,this%nm
     xa(n)=this%xa0+this%dxa*(n-1)
   enddo

   do i=1-this%ib,this%im+this%ib
     xf(i)=this%xf0+this%dxf*(i-1)
   enddo

   do m=1,this%mm
     ya(m)=this%ya0+this%dya*(m-1)
   enddo

   do j=1-this%jb,this%jm+this%jb
     yf(j)=this%yf0+this%dyf*(j-1)
   enddo

!
! Find iref and jref
!
   do n=1,this%nm
     do i=1-this%ib,this%im+this%ib-1
       if( xa(n)< xf(i)) then
         this%iref(n)=i-2
         exit
       endif
     enddo
   enddo

   do m=1,this%mm
     do j=1-this%jb,this%jm+this%jb-1
       if(ya(m) < yf(j)) then
         this%jref(m)=j-2
         exit
       endif
     enddo
   enddo

!ddreal(r_kind), dimension(1-this%ib:this%im+this%ib):: xf
write(6,*)"thinkdeb 0 ",1-this%ib, ' ',this%im+this%ib,this%nm

   do n=1,this%nm
     write(6,*)'thinkdeb n iref ',n,this%iref(n)
     i=this%iref(n)
     x1=xf(i)
     x2=xf(i+1)
     x3=xf(i+2)
     x4=xf(i+3)
     x = xa(n)
       x1x = x1-x   
       x2x = x2-x   
       x3x = x3-x   
       x4x = x4-x   
       rx2x1 = 1./(x2-x1)
       rx3x1 = 1./(x3-x1)
       rx4x1 = 1./(x4-x1)
       rx3x2 = 1./(x3-x2)
       rx4x2 = 1./(x4-x2)
       rx4x3 = 1./(x4-x3)
     CFL1 = x2x*x3x*rx2x1*rx3x1
     CFL2 =-x1x*x3x*rx2x1*rx3x2
     CFL3 = x1x*x2x*rx3x1*rx3x2
     CLL = x3x*rx3x2
     CFR1 = x3x*x4x*rx3x2*rx4x2
     CFR2 =-x2x*x4x*rx3x2*rx4x3
     CFR3 = x2x*x3x*rx4x2*rx4x3
     CRR =-x2x*rx3x2
       this%cx0(n)=CFL1*CLL
       this%cx1(n)=CFL2*CLL+CFR1*CRR
       this%cx2(n)=CFL3*CLL+CFR2*CRR
       this%cx3(n)=CFR3*CRR
   enddo

   do m=1,this%mm
     j=this%jref(m)
     y1=yf(j)
     y2=yf(j+1)
     y3=yf(j+2)
     y4=yf(j+3)
     y = ya(m)
       y1y = y1-y   
       y2y = y2-y   
       y3y = y3-y   
       y4y = y4-y   
       ry2y1 = 1./(y2-y1)
       ry3y1 = 1./(y3-y1)
       ry4y1 = 1./(y4-y1)
       ry3y2 = 1./(y3-y2)
       ry4y2 = 1./(y4-y2)
       ry4y3 = 1./(y4-y3)
     CFL1 = y2y*y3y*ry2y1*ry3y1
     CFL2 =-y1y*y3y*ry2y1*ry3y2
     CFL3 = y1y*y2y*ry3y1*ry3y2
     CLL = y3y*ry3y2
     CFR1 = y3y*y4y*ry3y2*ry4y2
     CFR2 =-y2y*y4y*ry3y2*ry4y3
     CFR3 = y2y*y3y*ry4y2*ry4y3
     CRR =-y2y*ry3y2
       this%cy0(m)=CFL1*CLL
       this%cy1(m)=CFL2*CLL+CFR1*CRR
       this%cy2(m)=CFL3*CLL+CFR2*CRR
       this%cy3(m)=CFR3*CRR
   enddo

 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_mg_coef

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     module                   subroutine lwq_vertical_coef &
!***********************************************************************
!                                                                      !
!  Prepare coeficients for vetical mapping between:                    !
!  analysis grid vertical resolution (nm) and                          !
!  generation one of filter grid vertical resoluition (im)             !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(this,nm_in,im_in,c1,c2,c3,c4,iref_out)
implicit none
class(mg_intstate_type),target::this

integer(i_kind), intent(in):: nm_in,im_in
real(r_kind), dimension(1:nm_in), intent(out):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(out):: iref_out

real(r_kind), dimension(1:nm_in):: y
real(r_kind), dimension(0:im_in+1):: x
real(r_kind):: dy,x1,x2,x3,x4,dx1,dx2,dx3,dx4 
real(r_kind):: dx13,dx23,dx24

integer(i_kind):: i,n
!-----------------------------------------------------------------------

   do i=0,im_in+1
     x(i)=(i-1)*1.
   enddo

    dy = 1.*(im_in-1)/(nm_in-1)
  do n=1,nm_in
    y(n)=(n-1)*dy
  enddo
    y(nm_in)=x(im_in)
 
  do n=2,nm_in-1
    i = y(n)+1
      x1 = x(i-1)
      x2 = x(i)
      x3 = x(i+1)
      x4 = x(i+2)
    iref_out(n)=i
      dx1 = y(n)-x1
      dx2 = y(n)-x2
      dx3 = y(n)-x3
      dx4 = y(n)-x4
      dx13 = dx1*dx3
      dx23 = 0.5*dx2*dx3
      dx24 = dx2*dx4
    c1(n) = -dx23*dx3
    c2(n) =  (    dx13+0.5*dx24)*dx3
    c3(n) = -(0.5*dx13+    dx24)*dx2
    c4(n) = dx23*dx2

    if(iref_out(n)==1) then
      c3(n)=c3(n)+c1(n)
      c1(n)=0.
    endif
    if(iref_out(n)==im_in-1) then
      c2(n)=c2(n)+c4(n)
      c4(n)=0.
    endif
  enddo
     iref_out(1)=1; c1(1)=0.; c2(1)=1.; c3(1)=0.; c4(1)=0.
     iref_out(nm_in)=im_in; c1(nm_in)=0.; c2(nm_in)=1.; c3(nm_in)=0.; c4(n)=0.
     
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_coef                            

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            module            subroutine lwq_vertical_adjoint                 &
!***********************************************************************
!                                                                      !
!  Direct linerly weighted quadratic adjoint interpolation in vertical !
!  from reslution nm to resolution im                                  !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(this,nm_in,km_in,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,w,f)
implicit none
!-----------------------------------------------------------------------
class(mg_intstate_type),target::this
integer(i_kind), intent(in):: nm_in,km_in,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm_in), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(in):: kref
real(r_kind), dimension(1:nm_in,imin:imax,jmin:jmax), intent(in):: w
real(r_kind), dimension(1:km_in,imin:imax,jmin:jmax), intent(out):: f
integer(i_kind):: k,n
!-----------------------------------------------------------------------
  f = 0.
do n=2,nm_in-1
  k = kref(n)
  if( k==1 ) then
    f(1,:,:) = f(1,:,:)+c2(n)*w(n,:,:)
    f(2,:,:) = f(2,:,:)+c3(n)*w(n,:,:)
    f(3,:,:) = f(3,:,:)+c4(n)*w(n,:,:)
  elseif &
    ( k==km_in-1) then
    f(km_in-2,:,:) = f(km_in-2,:,:)+c1(n)*w(n,:,:)
    f(km_in-1,:,:) = f(km_in-1,:,:)+c2(n)*w(n,:,:)
    f(km_in  ,:,:) = f(km_in  ,:,:)+c3(n)*w(n,:,:)
  elseif( k==km_in) then
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
  else
    f(k-1,:,:) = f(k-1,:,:)+c1(n)*w(n,:,:)
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
    f(k+1,:,:) = f(k+1,:,:)+c3(n)*w(n,:,:)
    f(k+2,:,:) = f(k+2,:,:)+c4(n)*w(n,:,:)
  endif
enddo
    f(1,:,:)=f(1,:,:)+w(1,:,:)
    f(km_in,:,:)=f(km_in,:,:)+w(nm_in,:,:)

!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_adjoint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 module       subroutine lwq_vertical_direct                  &
!***********************************************************************
!                                                                      !
!  Linerly weighted direct quadratic interpolation in vertical         !
!  from reslouion im to resolution nm                                  !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(this,km_in,nm_in,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,f,w)
implicit none
!-----------------------------------------------------------------------
class(mg_intstate_type),target::this
integer(i_kind), intent(in):: km_in,nm_in,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm_in), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(in):: kref
real(r_kind), dimension(1:km_in,imin:imax,jmin:jmax), intent(in):: f
real(r_kind), dimension(1:nm_in,imin:imax,jmin:jmax), intent(out):: w
integer(i_kind):: k,n
!-----------------------------------------------------------------------
do n=2,nm_in-1
  k = kref(n)
  if( k==1 ) then
    w(n,:,:) =             c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  elseif &
    ( k==km_in-1) then
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)
  elseif &
    ( k==km_in)   then
    w(n,:,:) =                 c2(n)*f(k,:,:)
  else
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,: )+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  endif
enddo
    w(1,:,:)=f(1,:,:)
    w(nm_in,:,:)=f(km_in,:,:)
    
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                      module  subroutine lsqr_direct_offset                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(km,i0-ib:im+ib,j0-jb:jm+jb) perform          !
! direct interpolations to get target array W(km,1:nm,1:mm)            !
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(this,V_in,W,km_in)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,1-this%ib:this%im+this%ib,1-this%jb:this%jm+this%jb), intent(in):: V_in
real(r_kind), dimension(km_in,1:this%nm,1:this%mm),intent(out):: W  

real(r_kind), dimension(km_in,1:this%nm,this%j0-this%jb:this%jm+this%jb):: VX
integer(i_kind):: i,j,n,m
real(r_kind),dimension(km_in):: v0,v1,v2,v3     
!-----------------------------------------------------------------------


   do j=this%j0-this%jb,this%jm+this%jb
   do n=1,this%nm
       i = this%iref(n)
     v0(:)=V_in(:,i  ,j)
     v1(:)=V_in(:,i+1,j)
     v2(:)=V_in(:,i+2,j)
     v3(:)=V_in(:,i+3,j)
     VX(:,n,j) = this%cx0(n)*v0(:)+this%cx1(n)*v1(:)+this%cx2(n)*v2(:)+this%cx3(n)*v3(:)
   enddo
   enddo

   do m=1,this%mm
     j = this%jref(m)
   do n=1,this%nm
     v0(:)=VX(:,n,j  ) 
     v1(:)=VX(:,n,j+1) 
     v2(:)=VX(:,n,j+2) 
     v3(:)=VX(:,n,j+3) 
     W(:,n,m) =  this%cy0(m)*v0(:)+this%cy1(m)*v1(:)+this%cy2(m)*v2(:)+this%cy3(m)*v3(:)
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_direct_offset

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           module             subroutine lsqr_adjoint_offset                  &
!***********************************************************************
!                                                                      !
! Given a target array W(km,0:nm,0:mm) perform adjoint                 !
! interpolations to get source array V(km,i0-ib:im+ib,j0-jb:jm+jb)     !
! using two passes of 1d interpolator                                  !
!                      - offset version -                              !
!                                                                      !
!***********************************************************************
(this,W,V_out,km_in)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
integer(i_kind):: km_in
real(r_kind), dimension(km_in,1:this%nm,1:this%mm),intent(in):: W  
real(r_kind), dimension(km_in,1-this%ib:this%im+this%ib,1-this%jb:this%jm+this%jb), intent(out):: V_out
real(r_kind), dimension(km_in,1:this%nm,this%j0-this%jb:this%jm+this%jb):: VX
integer(i_kind):: i,j,n,m,l,k
integer(i_kind):: ip1,ip2,ip3
integer(i_kind):: jp1,jp2,jp3
!-----------------------------------------------------------------------

   V_out(:,:,:) = 0.

   VX(:,:,:)=0.

   do m=1,this%mm
     j = this%jref(m)
     jp1=j+1
     jp2=j+2
     jp3=j+3
   do n=1,this%nm
     VX(:,n,j  ) = VX(:,n,j  )+W(:,n,m)*this%cy0(m)
     VX(:,n,jp1) = VX(:,n,jp1)+W(:,n,m)*this%cy1(m)
     VX(:,n,jp2) = VX(:,n,jp2)+W(:,n,m)*this%cy2(m)
     VX(:,n,jp3) = VX(:,n,jp3)+W(:,n,m)*this%cy3(m)
   enddo
   enddo
 

   do j=this%j0-this%jb,this%jm+this%jb
   do n=1,this%nm
     i = this%iref(n)
     ip1=i+1
     ip2=i+2
     ip3=i+3

     V_out(:,i  ,j) = V_out(:,i  ,j)+VX(:,n,j)*this%cx0(n)
     V_out(:,ip1,j) = V_out(:,ip1,j)+VX(:,n,j)*this%cx1(n)
     V_out(:,ip2,j) = V_out(:,ip2,j)+VX(:,n,j)*this%cx2(n)
     V_out(:,ip3,j) = V_out(:,ip3,j)+VX(:,n,j)*this%cx3(n)
   enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_offset

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule  mg_interpolate
