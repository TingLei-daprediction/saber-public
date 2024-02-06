!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_generations
!***********************************************************************
!                                                                      !
!  Contains procedures that include differrent generations             !
!                       - offset version -
!                                                                      !
!                                                     M. Rancic (2022) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter, only: i0,j0,im,jm,imL,jmL,hx,hy,gm
use mg_parameter, only: km,kmh,kmf,Fimax,Fjmax,FimaxL,FjmaxL
!use mpimod, only: mype     ! << for GSI  >>
use mg_mppstuff, only: mype
use mg_mppstuff, only: my_hgen,l_hgen,barrierMPI,finishMPI,Fimax,Fjmax
use mg_bocos, only: boco_2d,bocoT_2d
use mg_bocos, only: upsend_all,downsend_all
use mg_intstate, only: a_diff_h,b_diff_h
use mg_intstate, only: a_diff_f,b_diff_f
use mg_intstate, only: p_coef,q_coef
use mg_intstate, only: a_coef,b_coef
!TEST
use, intrinsic:: ieee_arithmetic
!TEST

public upsending_all
public downsending_all

public upsending2_all
public downsending2_all

public differencing_all

private adjoint_all
private direct_all

private adjoint2_all
private direct2_all


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsending_all                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(in):: V
real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(out):: H

real(r_kind),dimension(km,i0-2:imL+2,j0-2:jmL+2):: V_INT
real(r_kind),dimension(km,i0-2:imL+2,j0-2:jmL+2):: H_INT
integer(i_kind):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

        call adjoint_all(V(1:km,1:im,1:jm),V_INT,km,1) 

        call bocoT_2d(V_INT,km,imL,jmL,2,2)

        call upsend_all(V_INT(1:km,1:imL,1:jmL),H,km)
!
! From generation 2 sequentially to higher generations
!
  do g=2,gm-1 

    if(g==my_hgen) then
        call adjoint_all(H(1:km,1:im,1:jm),H_INT,km,g) 
    endif

        call bocoT_2d(H_INT,km,imL,jmL,2,2,FimaxL,FjmaxL,g,g)

        call upsend_all(H_INT(1:km,1:imL,1:jmL),H,km,g,g+1)

  end do    


!-----------------------------------------------------------------------
                        endsubroutine upsending_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsending_all                      &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(H,V)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(km,i0-hx:im+hx,j0-hy:jm+hy),intent(inout):: H
real(r_kind),dimension(km,i0-hx:im+hx,j0-hy:jm+hy),intent(inout):: V
real(r_kind),dimension(km,i0-2:imL+2,j0-2:jmL+2):: H_INT
real(r_kind),dimension(km,i0-2:imL+2,j0-2:jmL+2):: V_INT
real(r_kind),dimension(km,i0:im,j0:jm):: H_PROX
real(r_kind),dimension(km,i0:im,j0:jm):: V_PROX
integer(i_kind):: g,l,k
integer(i_kind):: iL,jL,i,j
!-----------------------------------------------------------------------
!
! Upper generations
!
    do g=gm,3,-1

        call downsend_all(H(1:km,i0:im,j0:jm),H_INT(1:km,1:imL,1:jmL),km,g,g-1)
        call boco_2d(H_INT,km,imL,jmL,2,2,FimaxL,FjmaxL,g-1,g-1)

      if(my_hgen==g-1) then
        call direct_all(H_INT,H_PROX,km,g-1)
        H(1:km,1:im,1:jm)=H     (1:km,i0:im,j0:jm)                    &
                         +H_PROX(1:km,i0:im,j0:jm)
      endif

    enddo

!
! From geneartion 2 to generation 1
!

        call downsend_all(H(1:km,i0:im,j0:jm),V_INT(1:km,1:imL,1:jmL),km)
          H(:,:,:)=0.

        call boco_2d(V_INT,km,imL,jmL,2,2)

        call direct_all(V_INT,V_PROX,km,1)

          V(1:km,i0:im,j0:jm)=V     (1:km,i0:im,j0:jm)                 &
                             +V_PROX(1:km,i0:im,j0:jm)

!-----------------------------------------------------------------------
                        endsubroutine downsending_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsending2_all                       &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(in):: V
real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(out):: H

real(r_kind),dimension(km,0:imL+1,0:jmL+1):: V_INT
real(r_kind),dimension(km,0:imL+1,0:jmL+1):: H_INT
integer(i_kind):: g,L
!-----------------------------------------------------------------------
!
! From generation 1 to generation 2
!

        call adjoint2_all(V(1:km,1:im,1:jm),V_INT,km,1) 

        call bocoT_2d(V_INT,km,imL,jmL,1,1)

        call upsend_all(V_INT(1:km,1:imL,1:jmL),H,km)
!
! From generation 2 sequentially to higher generations
!
  do g=2,gm-1 

    if(g==my_hgen) then
        call adjoint2_all(H(1:km,1:im,1:jm),H_INT,km,g) 
    endif

        call bocoT_2d(H_INT,km,imL,jmL,1,1,FimaxL,FjmaxL,g,g)

        call upsend_all(H_INT(1:km,1:imL,1:jmL),H,km,g,g+1)

  end do    


!-----------------------------------------------------------------------
                        endsubroutine upsending2_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsending2_all                     &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(H,V)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(inout):: H
real(r_kind),dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(inout):: V
real(r_kind),dimension(km,0:imL+1,0:jmL+1):: H_INT
real(r_kind),dimension(km,0:imL+1,0:jmL+1):: V_INT
real(r_kind),dimension(km,1:im,1:jm):: H_PROX
real(r_kind),dimension(km,1:im,1:jm):: V_PROX
integer(i_kind):: g,l,k
integer(i_kind):: iL,jL,i,j
!-----------------------------------------------------------------------
!
! Upper generations
!
    do g=gm,3,-1

        call downsend_all(H(1:km,1:im,1:jm),H_INT(1:km,1:imL,1:jmL),km,g,g-1)
        call boco_2d(H_INT,km,imL,jmL,1,1,FimaxL,FjmaxL,g-1,g-1)

      if(my_hgen==g-1) then
        call direct2_all(H_INT,H_PROX,km,g-1)
        H(1:km,1:im,1:jm)=H     (1:km,1:im,1:jm)                        &
                         +H_PROX(1:km,1:im,1:jm)
      endif

    enddo

!
! From generation 2 to generation 1
!

        call downsend_all(H(1:km,1:im,1:jm),V_INT(1:km,1:imL,1:jmL),km)
          H(:,:,:)=0.

        call boco_2d(V_INT,km,imL,jmL,1,1)

        call direct2_all(V_INT,V_PROX,km,1)

          V(1:km,1:im,1:jm)=V     (1:km,1:im,1:jm)                      &
                           +V_PROX(1:km,1:im,1:jm)

!-----------------------------------------------------------------------
                        endsubroutine downsending2_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine differencing_all                     &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(V,H)
!-----------------------------------------------------------------------
implicit none

real(r_kind),dimension(kmf,i0-hx:im+hx,j0-hy:jm+hy),intent(inout):: V
real(r_kind),dimension(kmh,i0-hx:im+hx,j0-hy:jm+hy),intent(inout):: H
real(r_kind),dimension(kmf,i0-1:im, j0  :jm):: DIFX
real(r_kind),dimension(kmf,i0  :im ,j0-1:jm):: DIFY
real(r_kind),dimension(kmh,i0-1:im, 0   :jm):: DIFXH
real(r_kind),dimension(kmh,i0  :im ,j0-1:jm):: DIFYH
integer(i_kind):: i,j,l,k,imx,jmx
!-----------------------------------------------------------------------

     do j=j0,jm
     do i=i0-1,im
       DIFX(:,i,j)=V(:,i+1,j)-V(:,i,j)
     enddo
     enddo
     do j=j0-1,jm
     do i=i0,im
       DIFY(:,i,j)=V(:,i,j+1)-V(:,i,j)
     enddo
     enddo


     do j=j0,jm
     do i=i0,im
       V(:,i,j)=a_diff_f(:,i,j)*V(:,i,j)                      &
               -b_diff_f(:,i,j)*(DIFX(:,i,j)-DIFX(:,i-1,j)    &
                                +DIFY(:,i,j)-DIFY(:,i,j-1))   
     enddo
     enddo

if(l_hgen) then

!  imx = Fimax(my_hgen)
!  jmx = Fjmax(my_hgen)

   imx = im
   jmx = jm

     do j=j0,jmx
     do i=i0-1,imx
       DIFXH(:,i,j)=H(:,i+1,j)-H(:,i,j)
     enddo
     enddo
     do j=j0-1,jmx
     do i=i0,imx
       DIFYH(:,i,j)=H(:,i,j+1)-H(:,i,j)
     enddo
     enddo

     do j=j0,jmx
     do i=i0,imx
        H(:,i,j)=a_diff_h(:,i,j)*H(:,i,j)                          &
                -b_diff_h(:,i,j)*(DIFXH(:,i,j)-DIFXH(:,i-1,j)      &
                                 +DIFYH(:,i,j)-DIFYH(:,i,j-1))  
     enddo
     enddo

endif

!-----------------------------------------------------------------------
                        endsubroutine differencing_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint_all                          &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           ! 
!                                                                      !
!***********************************************************************
(F,W,km,g)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: g 
integer(i_kind),intent(in):: km
real(r_kind), dimension(km,i0:im,j0:jm), intent(in):: F
real(r_kind), dimension(km,i0-2:imL+2,j0-2:jmL+2), intent(out):: W
real(r_kind), dimension(km,i0:im,j0-2:jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 3)
!
     W_AUX(:,:,:)= 0.

  do j=jm,2,-2
    jL = j/2
    do i=im,1,-1
      W_AUX(:,i,jL+2)=W_AUX(:,i,jL+2)+p_coef(4)*F(:,i,j)
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+p_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+p_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+p_coef(1)*F(:,i,j)
    enddo
  enddo
!
! 2)
!
  do j=jm-1,1,-2
    jL=j/2
    do i=im,1,-1
      W_AUX(:,i,jL+2)=W_AUX(:,i,jL+2)+q_coef(4)*F(:,i,j)
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+q_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+q_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+q_coef(1)*F(:,i,j)
    enddo
  enddo

    W(:,:,:)=0.
!
! 1)
!
  do jL=jmL+2,-1,-1
    do i=im-1,1,-2
    iL = i/2
      W(:,iL+2,jL)=W(:,iL+2,jL)+q_coef(4)*W_AUX(:,i,jL)
      W(:,iL+1,jL)=W(:,iL+1,jL)+q_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+q_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+q_coef(1)*W_AUX(:,i,jL)
    enddo
    do i=im,2,-2
    iL=i/2
      W(:,iL+2,jL)=W(:,iL+2,jL)+p_coef(4)*W_AUX(:,i,jL)
      W(:,iL+1,jL)=W(:,iL+1,jL)+p_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+p_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+p_coef(1)*W_AUX(:,i,jL)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine adjoint_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine direct_all                           &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(W,F,km,g)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km
real(r_kind), dimension(km,i0-2:imL+2,j0-2:jmL+2), intent(in):: W
real(r_kind), dimension(km,i0:im,j0:jm), intent(out):: F
real(r_kind), dimension(km,i0:im,j0-2:jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------

!
! 1)
!
   do jL=-1,jmL+2
     do i=1,im-1,2
       iL=i/2
         W_AUX(:,i,jL)=q_coef(1)*W(:,iL-1,jL)+q_coef(2)*W(:,iL  ,jL)              &
                      +q_coef(3)*W(:,iL+1,jL)+q_coef(4)*W(:,iL+2,jL)
     enddo
     do i=2,im,2
       iL=i/2
         W_AUX(:,i,jL)=p_coef(1)*W(:,iL-1,jL)+p_coef(2)*w(:,iL  ,jL)              &
                      +p_coef(3)*W(:,iL+1,jL)+p_coef(4)*W(:,iL+2,jL)
     enddo
   enddo
!
! 2)
!
   do j=1,jm-1,2
     jL=j/2
     do i=1,im
       F(:,i,j)=q_coef(1)*W_AUX(:,i,jL-1)+q_coef(2)*W_AUX(:,i,jL  )               &
               +q_coef(3)*W_AUX(:,i,jL+1)+q_coef(4)*W_AUX(:,i,jL+2)
     enddo
   enddo
!
! 3)
!
   do j=2,jm,2
     jL=j/2
     do i=1,im
       F(:,i,j)=p_coef(1)*W_AUX(:,i,jL-1)+p_coef(2)*W_AUX(:,i,jL  )               &
               +p_coef(3)*W_AUX(:,i,jL+1)+p_coef(4)*W_AUX(:,i,jL+2)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine direct_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine adjoint2_all                         &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using quadratics interpolations                                    !
!                         - offset version -                           ! 
!                                                                      !
!***********************************************************************
(F,W,km,g)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: g 
integer(i_kind),intent(in):: km
real(r_kind), dimension(km,1:im,1:jm), intent(in):: F
real(r_kind), dimension(km,0:imL+1,0:jmL+1), intent(out):: W
real(r_kind), dimension(km,1:im,0:jmL+2):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 3)
!
     W_AUX(:,:,:)= 0.

  do j=jm,2,-2
    jL = j/2
    do i=im,1,-1
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+b_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+b_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+b_coef(1)*F(:,i,j)
    enddo
  enddo
!
! 2)
!
  do j=jm-1,1,-2
    jL=(j+1)/2
    do i=im,1,-1
      W_AUX(:,i,jL+1)=W_AUX(:,i,jL+1)+a_coef(3)*F(:,i,j)
      W_AUX(:,i,jL  )=W_AUX(:,i,jL  )+a_coef(2)*F(:,i,j)
      W_AUX(:,i,jL-1)=W_AUX(:,i,jL-1)+a_coef(1)*F(:,i,j)
    enddo
  enddo

    W(:,:,:)=0.
!
! 1)
!
  do jL=jmL+1,0,-1
    do i=im-1,1,-2
    iL = (i+1)/2
      W(:,iL+1,jL)=W(:,iL+1,jL)+a_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+a_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+a_coef(1)*W_AUX(:,i,jL)
    enddo
    do i=im,2,-2
    iL=i/2
      W(:,iL+1,jL)=W(:,iL+1,jL)+b_coef(3)*W_AUX(:,i,jL)
      W(:,iL  ,jL)=W(:,iL  ,jL)+b_coef(2)*W_AUX(:,i,jL)
      W(:,iL-1,jL)=W(:,iL-1,jL)+b_coef(1)*W_AUX(:,i,jL)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine adjoint2_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine direct2_all                          &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using quadratic interpolations                                     !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(W,F,km,g)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km
real(r_kind), dimension(km,0:imL+1,0:jmL+1), intent(in):: W
real(r_kind), dimension(km,1:im,1:jm), intent(out):: F
real(r_kind), dimension(km,1:im,0:jmL+1):: W_AUX
integer(i_kind):: i,j,iL,jL
!-----------------------------------------------------------------------
!
! 1)
!
   do jL=0,jmL+1
     do i=1,im-1,2
       iL=(i+1)/2
         W_AUX(:,i,jL)=a_coef(1)*W(:,iL-1,jL)+a_coef(2)*W(:,iL  ,jL)    &
                      +a_coef(3)*W(:,iL+1,jL)
     enddo
     do i=2,im,2
       iL=i/2
         W_AUX(:,i,jL)=b_coef(1)*W(:,iL-1,jL)+b_coef(2)*w(:,iL  ,jL)    &
                      +b_coef(3)*W(:,iL+1,jL)
     enddo
   enddo
!
! 2)
!
   do j=1,jm-1,2
     jL=(j+1)/2
     do i=1,im
       F(:,i,j)=a_coef(1)*W_AUX(:,i,jL-1)+a_coef(2)*W_AUX(:,i,jL  )     &
               +a_coef(3)*W_AUX(:,i,jL+1)
     enddo
   enddo
!
! 3)
!
   do j=2,jm,2
     jL=j/2
     do i=1,im
       F(:,i,j)=b_coef(1)*W_AUX(:,i,jL-1)+b_coef(2)*W_AUX(:,i,jL  )     &
               +b_coef(3)*W_AUX(:,i,jL+1)
     enddo
   enddo

!-----------------------------------------------------------------------
                        endsubroutine direct2_all

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_generations
