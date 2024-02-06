!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        submodule(mg_intstate)  mg_generations
!***********************************************************************
!                                                                      !
!  Contains procedures that include differrent generations             !
!                       - offset version -
!                                                                      !
!                                                     M. Rancic (2022) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_timers
!TEST
use, intrinsic:: ieee_arithmetic
!TEST


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
               module         subroutine upsending_all                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!                                                                      !
!***********************************************************************
(this,V,H,lquart)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
logical, intent(in):: lquart
!-----------------------------------------------------------------------

        if(lquart) then
           call this%upsending2(V,H) 
        else
           call this%upsending(V,H) 
        endif
         

!-----------------------------------------------------------------------
                        endsubroutine upsending_all 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         module               subroutine downsending_all                      &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(this,H,V,lquart)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V

logical, intent(in):: lquart
!-----------------------------------------------------------------------

        if(lquart) then
           call this%downsending2(H,V) 
        else
           call this%downsending(H,V) 
        endif

!-----------------------------------------------------------------------
                        endsubroutine downsending_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           module             subroutine weighting_all                        &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H,lhelm)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H

logical, intent(in):: lhelm
!-----------------------------------------------------------------------

        if(lhelm) then
           call this%weighting_helm(V,H) 
        else
           call this%weighting(V,H) 
        endif

!-----------------------------------------------------------------------
                        endsubroutine weighting_all 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          module              subroutine upsending                            &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H

real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: V_INT
real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: H_INT
integer(i_kind):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

        call this%adjoint(V(1:this%km,1:this%im,1:this%jm),V_INT,this%km,1) 

        call this%bocoT_2d(V_INT,this%km,this%imL,this%jmL,2,2)

        call this%upsend_all(V_INT(1:this%km,1:this%imL,1:this%jmL),H,this%km)
!
! From generation 2 sequentially to higher generations
!
  do g=2,this%gm-1 

    if(g==this%my_hgen) then
        call this%adjoint(H(1:this%km,1:this%im,1:this%jm),H_INT,this%km,g) 
    endif

        call this%bocoT_2d(H_INT,this%km,this%imL,this%jmL,2,2,this%FimaxL,this%FjmaxL,g,g)

        call this%upsend_all(H_INT(1:this%km,1:this%imL,1:this%jmL),H,this%km,g,g+1)

  end do    


!-----------------------------------------------------------------------
                        endsubroutine upsending

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
           module             subroutine downsending                          &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(this,H,V)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V

real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: H_INT
real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: V_INT
real(r_kind),dimension(this%km,this%i0:this%im,this%j0:this%jm):: H_PROX
real(r_kind),dimension(this%km,this%i0:this%im,this%j0:this%jm):: V_PROX
integer(i_kind):: g,l,k
integer(i_kind):: iL,jL,i,j
!-----------------------------------------------------------------------
!
! Upper generations
!
    do g=this%gm,3,-1

        call this%downsend_all(H(1:this%km,this%i0:this%im,this%j0:this%jm),H_INT(1:this%km,1:this%imL,1:this%jmL),this%km,g,g-1)
        call this%boco_2d(H_INT,this%km,this%imL,this%jmL,2,2,this%FimaxL,this%FjmaxL,g-1,g-1)

      if(this%my_hgen==g-1) then
        call this%direct1(H_INT,H_PROX,this%km,g-1)
        H(1:this%km,1:this%im,1:this%jm)=H     (1:this%km,this%i0:this%im,this%j0:this%jm)                    &
                         +H_PROX(1:this%km,this%i0:this%im,this%j0:this%jm)
      endif

    enddo

!
! From geneartion 2 to generation 1
!

        call this%downsend_all(H(1:this%km,this%i0:this%im,this%j0:this%jm),V_INT(1:this%km,1:this%imL,1:this%jmL),this%km)
          H(:,:,:)=0.

        call this%boco_2d(V_INT,this%km,this%imL,this%jmL,2,2)

        call this%direct1(V_INT,V_PROX,this%km,1)

          V(1:this%km,this%i0:this%im,this%j0:this%jm)=V     (1:this%km,this%i0:this%im,this%j0:this%jm)                 &
                             +V_PROX(1:this%km,this%i0:this%im,this%j0:this%jm)

!-----------------------------------------------------------------------
                        endsubroutine downsending

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                 module       subroutine upsending2                           &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H


real(r_kind),dimension(this%km,0:this%imL+1,0:this%jmL+1):: V_INT
real(r_kind),dimension(this%km,0:this%imL+1,0:this%jmL+1):: H_INT
integer(i_kind):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

        call this%adjoint2(V(1:this%km,1:this%im,1:this%jm),V_INT,this%km,1) 

        call this%bocoT_2d(V_INT,this%km,this%imL,this%jmL,1,1)

        call this%upsend_all(V_INT(1:this%km,1:this%imL,1:this%jmL),H,this%km)
!
! From generation 2 sequentially to higher generations
!
  do g=2,this%gm-1 

    if(g==this%my_hgen) then
        call this%adjoint2(H(1:this%km,1:this%im,1:this%jm),H_INT,this%km,g) 
    endif

        call this%bocoT_2d(H_INT,this%km,this%imL,this%jmL,1,1,this%FimaxL,this%FjmaxL,g,g)

        call this%upsend_all(H_INT(1:this%km,1:this%imL,1:this%jmL),H,this%km,g,g+1)

  end do    


!-----------------------------------------------------------------------
                        endsubroutine upsending2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                    module    subroutine downsending2                         &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(this,H,V)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,0:this%imL+1,0:this%jmL+1):: H_INT
real(r_kind),dimension(this%km,0:this%imL+1,0:this%jmL+1):: V_INT
real(r_kind),dimension(this%km,1:this%im,1:this%jm):: H_PROX
real(r_kind),dimension(this%km,1:this%im,1:this%jm):: V_PROX
integer(i_kind):: g,l,k
integer(i_kind):: iL,jL,i,j
!-----------------------------------------------------------------------
!
! Upper generations
!
    do g=this%gm,3,-1

        call this%downsend_all(H(1:this%km,1:this%im,1:this%jm),H_INT(1:this%km,1:this%imL,1:this%jmL),this%km,g,g-1)
        call this%boco_2d(H_INT,this%km,this%imL,this%jmL,1,1,this%FimaxL,this%FjmaxL,g-1,g-1)

      if(this%my_hgen==g-1) then
        call this%direct2(H_INT,H_PROX,this%km,g-1)
        H(1:this%km,1:this%im,1:this%jm)=H     (1:this%km,1:this%im,1:this%jm)                        &
                         +H_PROX(1:this%km,1:this%im,1:this%jm)
      endif

    enddo

!
! From generation 2 to generation 1
!

        call this%downsend_all(H(1:this%km,1:this%im,1:this%jm),V_INT(1:this%km,1:this%imL,1:this%jmL),this%km)
          H(:,:,:)=0.

        call this%boco_2d(V_INT,this%km,this%imL,this%jmL,1,1)

        call this%direct2(V_INT,V_PROX,this%km,1)

          V(1:this%km,1:this%im,1:this%jm)=V     (1:this%km,1:this%im,1:this%jm)                      &
                           +V_PROX(1:this%km,1:this%im,1:this%jm)

!-----------------------------------------------------------------------
                        endsubroutine downsending2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                     module   subroutine weighting_helm                       &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,this%i0-1:this%im, this%j0  :this%jm):: DIFX
real(r_kind),dimension(this%km,this%i0  :this%im ,this%j0-1:this%jm):: DIFY
real(r_kind),dimension(this%km,this%i0-1:this%im, this%j0  :this%jm):: DIFXH
real(r_kind),dimension(this%km,this%i0  :this%im ,this%j0-1:this%jm):: DIFYH
integer(i_kind):: i,j,l,k,imx,jmx
!-----------------------------------------------------------------------

     do j=this%j0,this%jm
     do i=this%i0-1,this%im
       DIFX(:,i,j)=V(:,i+1,j)-V(:,i,j)
     enddo
     enddo
     do j=this%j0-1,this%jm
     do i=this%i0,this%im
       DIFY(:,i,j)=V(:,i,j+1)-V(:,i,j)
     enddo
     enddo


     do j=this%j0,this%jm
     do i=this%i0,this%im
       V(:,i,j)=this%a_diff_f(:,i,j)*V(:,i,j)                      &
               -this%b_diff_f(:,i,j)*(DIFX(:,i,j)-DIFX(:,i-1,j)    &
                                +DIFY(:,i,j)-DIFY(:,i,j-1))   
     enddo
     enddo

if(this%l_hgen) then

!  imx = Fimax(my_hgen)
!  jmx = Fjmax(my_hgen)

   imx = this%im
   jmx = this%jm

     do j=this%j0,jmx
     do i=this%i0-1,imx
       DIFXH(:,i,j)=H(:,i+1,j)-H(:,i,j)
     enddo
     enddo
     do j=this%j0-1,jmx
     do i=this%i0,imx
       DIFYH(:,i,j)=H(:,i,j+1)-H(:,i,j)
     enddo
     enddo

     do j=this%j0,jmx
     do i=this%i0,imx
        H(:,i,j)=this%a_diff_h(:,i,j)*H(:,i,j)                          &
                -this%b_diff_h(:,i,j)*(DIFXH(:,i,j)-DIFXH(:,i-1,j)      &
                                 +DIFYH(:,i,j)-DIFYH(:,i,j-1))  
     enddo
     enddo

endif

!-----------------------------------------------------------------------
                        endsubroutine weighting_helm

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                     module   subroutine weighting                            &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H

integer(i_kind):: i,j,l,k,imx,jmx
!-----------------------------------------------------------------------

     do j=this%j0,this%jm
     do i=this%i0,this%im
       V(:,i,j)=this%a_diff_f(:,i,j)*V(:,i,j)                      
     enddo
     enddo

if(this%l_hgen) then

   imx = this%im
   jmx = this%jm

     do j=this%j0,jmx
     do i=this%i0,imx
        H(:,i,j)=this%a_diff_h(:,i,j)*H(:,i,j)                          
     enddo
     enddo

endif


!-----------------------------------------------------------------------
                        endsubroutine weighting 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                     module   subroutine adjoint                              &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           ! 
!                                                                      !
!***********************************************************************
(this,F,W,km_in,g)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g 
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,this%i0:this%im,this%j0:this%jm), intent(in):: F
real(r_kind), dimension(km_in,this%i0-2:this%imL+2,this%j0-2:this%jmL+2), intent(out):: W
real(r_kind), dimension(km_in,this%i0:this%im,this%j0-2:this%jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 3)
!
     W_AUX(:,:,:)= 0.

  do j=this%jm,2,-2
    jL = j/2
    do i=this%im,1,-1
      W_AUX(:,i,jL+2)=W_AUX(:,i,jL+2)+this%p_coef(4)*F(:,i,j)
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+this%p_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+this%p_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+this%p_coef(1)*F(:,i,j)
    enddo
  enddo
!
! 2)
!
  do j=this%jm-1,1,-2
    jL=j/2
    do i=this%im,1,-1
      W_AUX(:,i,jL+2)=W_AUX(:,i,jL+2)+this%q_coef(4)*F(:,i,j)
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+this%q_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+this%q_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+this%q_coef(1)*F(:,i,j)
    enddo
  enddo

    W(:,:,:)=0.
!
! 1)
!
  do jL=this%jmL+2,-1,-1
    do i=this%im-1,1,-2
    iL = i/2
      W(:,iL+2,jL)=W(:,iL+2,jL)+this%q_coef(4)*W_AUX(:,i,jL)
      W(:,iL+1,jL)=W(:,iL+1,jL)+this%q_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+this%q_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+this%q_coef(1)*W_AUX(:,i,jL)
    enddo
    do i=this%im,2,-2
    iL=i/2
      W(:,iL+2,jL)=W(:,iL+2,jL)+this%p_coef(4)*W_AUX(:,i,jL)
      W(:,iL+1,jL)=W(:,iL+1,jL)+this%p_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+this%p_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+this%p_coef(1)*W_AUX(:,i,jL)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine adjoint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                    module    subroutine direct1                              &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,W,F,km_in,g)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,this%i0-2:this%imL+2,this%j0-2:this%jmL+2), intent(in):: W
real(r_kind), dimension(km_in,this%i0:this%im,this%j0:this%jm), intent(out):: F
real(r_kind), dimension(km_in,this%i0:this%im,this%j0-2:this%jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------

!
! 1)
!
   do jL=-1,this%jmL+2
     do i=1,this%im-1,2
       iL=i/2
         W_AUX(:,i,jL)=this%q_coef(1)*W(:,iL-1,jL)+this%q_coef(2)*W(:,iL  ,jL)              &
                      +this%q_coef(3)*W(:,iL+1,jL)+this%q_coef(4)*W(:,iL+2,jL)
     enddo
     do i=2,this%im,2
       iL=i/2
         W_AUX(:,i,jL)=this%p_coef(1)*W(:,iL-1,jL)+this%p_coef(2)*w(:,iL  ,jL)              &
                      +this%p_coef(3)*W(:,iL+1,jL)+this%p_coef(4)*W(:,iL+2,jL)
     enddo
   enddo
!
! 2)
!
   do j=1,this%jm-1,2
     jL=j/2
     do i=1,this%im
       F(:,i,j)=this%q_coef(1)*W_AUX(:,i,jL-1)+this%q_coef(2)*W_AUX(:,i,jL  )               &
               +this%q_coef(3)*W_AUX(:,i,jL+1)+this%q_coef(4)*W_AUX(:,i,jL+2)
     enddo
   enddo
!
! 3)
!
   do j=2,this%jm,2
     jL=j/2
     do i=1,this%im
       F(:,i,j)=this%p_coef(1)*W_AUX(:,i,jL-1)+this%p_coef(2)*W_AUX(:,i,jL  )               &
               +this%p_coef(3)*W_AUX(:,i,jL+1)+this%p_coef(4)*W_AUX(:,i,jL+2)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine direct1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                    module    subroutine adjoint2                             &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using quadratics interpolations                                    !
!                         - offset version -                           ! 
!                                                                      !
!***********************************************************************
(this,F,W,km_in,g)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g 
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm), intent(in):: F
real(r_kind), dimension(km_in,0:this%imL+1,0:this%jmL+1), intent(out):: W
real(r_kind), dimension(km_in,1:this%im,0:this%jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 3)
!
     W_AUX(:,:,:)= 0.

  do j=this%jm,2,-2
    jL = j/2
    do i=this%im,1,-1
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+this%b_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+this%b_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+this%b_coef(1)*F(:,i,j)
    enddo
  enddo
!
! 2)
!
  do j=this%jm-1,1,-2
    jL=(j+1)/2
    do i=this%im,1,-1
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+this%a_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+this%a_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+this%a_coef(1)*F(:,i,j)
    enddo
  enddo

    W(:,:,:)=0.
!
! 1)
!
  do jL=this%jmL+1,0,-1
    do i=this%im-1,1,-2
    iL = (i+1)/2
      W(:,iL+1,jL)=W(:,iL+1,jL)+this%a_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+this%a_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+this%a_coef(1)*W_AUX(:,i,jL)
    enddo
    do i=this%im,2,-2
    iL=i/2
      W(:,iL+1,jL)=W(:,iL+1,jL)+this%b_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+this%b_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+this%b_coef(1)*W_AUX(:,i,jL)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine adjoint2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                     module   subroutine direct2                              &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using quadratic interpolations                                     !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,W,F,km_in,g)
!-----------------------------------------------------------------------
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,0:this%imL+1,0:this%jmL+1), intent(in):: W
real(r_kind), dimension(km_in,1:this%im,1:this%jm), intent(out):: F
real(r_kind), dimension(km_in,1:this%im,0:this%jmL+1):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 1)
!
   do jL=0,this%jmL+1
     do i=1,this%im-1,2
       iL=(i+1)/2
         W_AUX(:,i,jL)=this%a_coef(1)*W(:,iL-1,jL)+this%a_coef(2)*W(:,iL  ,jL)    &
                      +this%a_coef(3)*W(:,iL+1,jL)
     enddo
     do i=2,this%im,2
       iL=i/2
         W_AUX(:,i,jL)=this%b_coef(1)*W(:,iL-1,jL)+this%b_coef(2)*w(:,iL  ,jL)    &
                      +this%b_coef(3)*W(:,iL+1,jL)
     enddo
   enddo
!
! 2)
!
   do j=1,this%jm-1,2
     jL=(j+1)/2
     do i=1,this%im
       F(:,i,j)=this%a_coef(1)*W_AUX(:,i,jL-1)+this%a_coef(2)*W_AUX(:,i,jL  )     &
               +this%a_coef(3)*W_AUX(:,i,jL+1)
     enddo
   enddo
!
! 3)
!
   do j=2,this%jm,2
     jL=j/2
     do i=1,this%im
       F(:,i,j)=this%b_coef(1)*W_AUX(:,i,jL-1)+this%b_coef(2)*W_AUX(:,i,jL  )     &
               +this%b_coef(3)*W_AUX(:,i,jL+1)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine direct2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end submodule mg_generations
