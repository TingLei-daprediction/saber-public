!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_intstate
!***********************************************************************
!                                                                      !
! Contains declarations and allocations of internal state variables    !
! use for filtering                                                    !
!                       - offset version -                             !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use jp_pkind2, only: fpi
!GSI use mpimod, only: mype
!use mg_entrymod, only: km2,km3,km
!GSI use berror, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
!clt use jp_pbfil,only: cholaspect
!use jp_pbfil,only: getlinesum
use jp_pbfil3, only: inimomtab,t22_to_3,tritform,t33_to_6,hextform
use mg_parameter,only: mg_parameter_type
!TEST
!use gridmod, only: lat1,lon1
!TEST
implicit none
type,extends( mg_parameter_type):: mg_intstate_type
real(r_kind), allocatable,dimension(:,:,:):: V
!
! Composite control variable on first generation o filter grid
!
real(r_kind), allocatable,dimension(:,:,:):: VALL
real(r_kind), allocatable,dimension(:,:,:):: HALL
!
! Composite control variable on high generations of filter grid
!
!
!FOR ADJOINT TEST
!
!real(r_kind), allocatable,dimension(:,:):: A
!real(r_kind), allocatable,dimension(:,:):: B
!real(r_kind), allocatable,dimension(:,:):: A0
!real(r_kind), allocatable,dimension(:,:):: B0
!
real(r_kind), allocatable,dimension(:,:,:):: a_diff_f
real(r_kind), allocatable,dimension(:,:,:):: a_diff_h
real(r_kind), allocatable,dimension(:,:,:):: b_diff_f
real(r_kind), allocatable,dimension(:,:,:):: b_diff_h

real(r_kind), allocatable,dimension(:,:):: p_eps
real(r_kind), allocatable,dimension(:,:):: p_del
real(r_kind), allocatable,dimension(:,:):: p_sig
real(r_kind), allocatable,dimension(:,:):: p_rho

real(r_kind), allocatable,dimension(:,:,:):: paspx
real(r_kind), allocatable,dimension(:,:,:):: paspy
real(r_kind), allocatable,dimension(:,:,:):: pasp1
real(r_kind), allocatable,dimension(:,:,:,:):: pasp2
real(r_kind), allocatable,dimension(:,:,:,:,:):: pasp3

real(r_kind), allocatable,dimension(:,:,:):: vpasp2
real(r_kind), allocatable,dimension(:,:,:):: hss2
real(r_kind), allocatable,dimension(:,:,:,:):: vpasp3
real(r_kind), allocatable,dimension(:,:,:,:):: hss3

real(r_kind), allocatable,dimension(:):: ssx
real(r_kind), allocatable,dimension(:):: ssy
real(r_kind), allocatable,dimension(:):: ss1
real(r_kind), allocatable,dimension(:,:):: ss2
real(r_kind), allocatable,dimension(:,:,:):: ss3

integer(fpi), allocatable,dimension(:,:,:):: dixs
integer(fpi), allocatable,dimension(:,:,:):: diys
integer(fpi), allocatable,dimension(:,:,:):: dizs

integer(fpi), allocatable,dimension(:,:,:,:):: dixs3
integer(fpi), allocatable,dimension(:,:,:,:):: diys3
integer(fpi), allocatable,dimension(:,:,:,:):: dizs3

integer(fpi), allocatable,dimension(:,:,:,:):: qcols

!real(r_kind), allocatable,dimension(:,:,:,:):: r_vol
!
!
! Composite stacked variable
!

!cltreal(r_kind), allocatable,dimension(:,:,:):: WORKA


integer(i_kind),allocatable,dimension(:):: iref,jref
integer(i_kind),allocatable,dimension(:):: Lref,Lref_h
real(r_kind),allocatable,dimension(:):: cvf1,cvf2,cvf3,cvf4
real(r_kind),allocatable,dimension(:):: cvh1,cvh2,cvh3,cvh4

real(r_kind),allocatable,dimension(:):: cx0,cx1,cx2,cx3
real(r_kind),allocatable,dimension(:):: cy0,cy1,cy2,cy3

real(r_kind),allocatable,dimension(:):: p_coef,q_coef
real(r_kind),allocatable,dimension(:):: a_coef,b_coef

real(r_kind),allocatable,dimension(:,:):: cf00,cf01,cf02,cf03           &
                                         ,cf10,cf11,cf12,cf13           &
                                         ,cf20,cf21,cf22,cf23           &
                                         ,cf30,cf31,cf32,cf33
!clt from interpolate.f90
contains
        procedure :: allocate_mg_intstate, def_mg_weights, init_mg_line
        procedure  ::def_offset_coef
        procedure  ::lsqr_mg_coef,lwq_vertical_coef 
        procedure  ::lwq_vertical_direct,lwq_vertical_adjoint , &
                 lsqr_direct_offset,  &
                 deallocate_mg_intstate, &
                 lsqr_adjoint_offset 
        generic ::boco_2d => boco_2d_g1,boco_2d_gh 
        generic ::boco_3d => boco_3d_g1,boco_3d_gh 
        generic ::bocoT_2d => bocoT_2d_g1,bocoT_2d_gh 
        generic ::bocoTx => bocoTx_2d_g1,bocoTx_2d_gh 
        generic ::bocoTy => bocoTy_2d_g1,bocoTy_2d_gh 
        generic ::bocoT_3d => bocoT_3d_g1,bocoT_3d_gh 
        generic ::bocox => bocox_2d_g1,bocox_2d_gh 
        generic ::bocoy => bocoy_2d_g1,bocoy_2d_gh 

        generic ::upsend_all=> upsend_all_g1 ,upsend_all_gh
        generic ::downsend_all=> downsend_all_g2 ,downsend_all_gh
        procedure:: upsend_all_g1 ,upsend_all_gh
        procedure:: downsend_all_g2 ,downsend_all_gh
        procedure:: boco_2d_g1,boco_2d_gh
        procedure:: boco_3d_g1,boco_3d_gh
        procedure :: bocoT_2d_g1,bocoT_2d_gh
        procedure :: bocoTx_2d_g1,bocoTx_2d_gh 
        procedure :: bocoTy_2d_g1,bocoTy_2d_gh
        procedure :: bocoT_3d_g1,bocoT_3d_gh 
        procedure ::  bocox_2d_g1,bocox_2d_gh 
        procedure  :: bocoy_2d_g1,bocoy_2d_gh 
!cltfrom mg_generation
        procedure:: upsending_all,downsending_all,weighting_all, &
               upsending,downsending,upsending2,downsending2, &
               weighting_helm,weighting ,adjoint,direct1, &
               adjoint2,direct2
!clt mg_filtering 
         procedure ::sup_vrbeta1T,sup_vrbeta1,sup_vrbeta3T,sup_vrbeta3

         procedure:: mg_filtering_rad1,mg_filtering_rad2,mg_filtering_rad3,&
                 mg_filtering_lin1,mg_filtering_lin2,mg_filtering_lin3, &
                 mg_filtering_fast
!clt from mg_transfer.f90
          procedure:: composite_to_stack,stack_to_composite
!clt from mg_entrymod
          procedure :: mg_initialize
          procedure ::mg_finalize
          procedure :: anal_to_filt_all,mg_filtering_procedure,filt_to_anal_all
end type mg_intstate_type
 interface 
 module subroutine  lsqr_mg_coef(this)
   import mg_intstate_type
 class(mg_intstate_type),target::this
 end subroutine
 module subroutine lwq_vertical_coef  &
(this,nm_in,im_in,c1,c2,c3,c4,iref_out)
   import mg_intstate_type
implicit none
 class(mg_intstate_type),target::this

integer(i_kind), intent(in):: nm_in,im_in
real(r_kind), dimension(1:nm_in), intent(out):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(out):: iref_out
 end subroutine

 module subroutine  lwq_vertical_direct &
(this,km_in,nm_in,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,f,w)
   import mg_intstate_type
implicit none
!-----------------------------------------------------------------------
class(mg_intstate_type),target::this
integer(i_kind), intent(in):: km_in,nm_in,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm_in), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(in):: kref
real(r_kind), dimension(1:km_in,imin:imax,jmin:jmax), intent(in):: f
real(r_kind), dimension(1:nm_in,imin:imax,jmin:jmax), intent(out):: w
 end subroutine
 module subroutine lwq_vertical_adjoint  &
(this,nm_in,km_in,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,w,f)
   import mg_intstate_type
implicit none
!-----------------------------------------------------------------------
class(mg_intstate_type),target::this
integer(i_kind), intent(in):: nm_in,km_in,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm_in), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm_in), intent(in):: kref
real(r_kind), dimension(1:nm_in,imin:imax,jmin:jmax), intent(in):: w
real(r_kind), dimension(1:km_in,imin:imax,jmin:jmax), intent(out):: f
 end subroutine
end interface 
interface 
 module subroutine def_offset_coef(this)
   import mg_intstate_type
 class(mg_intstate_type),target::this
 end subroutine
end interface 

interface 

         module    subroutine lsqr_direct_offset                   &
(this,V_in,W,km_in)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,1-this%ib:this%im+this%ib,1-this%jb:this%jm+this%jb), intent(in):: V_in
real(r_kind), dimension(km_in,1:this%nm,1:this%mm),intent(out):: W  

real(r_kind), dimension(km_in,1:this%nm,this%j0-this%jb:this%jm+this%jb):: VX
         end    subroutine lsqr_direct_offset                   

         module               subroutine lsqr_adjoint_offset                  &
(this,W,V_out,km_in)
!-----------------------------------------------------------------------
   import mg_intstate_type
   import i_kind,r_kind
implicit none
class(mg_intstate_type),target::this
integer(i_kind):: km_in
real(r_kind), dimension(km_in,1:this%nm,1:this%mm),intent(in):: W  
real(r_kind), dimension(km_in,1-this%ib:this%im+this%ib,1-this%jb:this%jm+this%jb), intent(out):: V_out
real(r_kind), dimension(km_in,1:this%nm,this%j0-this%jb:this%jm+this%jb):: VX
 end subroutine
!clt from mg_transfer.f90

 module  subroutine anal_to_filt_all(this,WORKA)
   import mg_intstate_type
   class(mg_intstate_type),target::this
   real (r_kind):: WORKA(this%km,this%n0:this%nm,this%m0:this%mm)

  end subroutine anal_to_filt_all
  module subroutine filt_to_anal_all (this,WORKA)
   import mg_intstate_type
    class(mg_intstate_type),target::this
   real (r_kind):: WORKA(this%km,this%n0:this%nm,this%m0:this%mm)

  end subroutine filt_to_anal_all 
  module subroutine stack_to_composite(this,ARR_ALL,A2D,A3D)
   import mg_intstate_type
   class(mg_intstate_type),target::this
real(r_kind),dimension(this%km ,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(in):: ARR_ALL
real(r_kind),dimension(this%km3,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy,this%lm),intent(out):: A3D
real(r_kind),dimension(this%km2,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy)   ,intent(out):: A2D


  end subroutine stack_to_composite 
module  subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(this,A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class(mg_intstate_type),target::this
real(r_kind),dimension(this%km2,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(in):: A2D
real(r_kind),dimension(this%km3,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy,this%lm),intent(in):: A3D
real(r_kind),dimension(this%km ,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),   intent(out):: ARR_ALL
  end subroutine composite_to_stack
  
 end interface
!clt for mg_bocos 
interface 
   module              subroutine boco_2d_g1                           &
(this,W,km_in,im_in,jm_in,nbx,nby)
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
   end               subroutine 
   module    subroutine boco_2d_gh                           &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
   end    subroutine boco_2d_gh 

             module           subroutine boco_3d_g1                           &
(this,W,km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in)
!-----------------------------------------------------------------------

implicit none

class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz
real(r_kind),dimension(km3_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz)    &
                      ,intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
             end            subroutine boco_3d_g1

            module      subroutine boco_3d_gh                           &
(this,W,km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------

implicit none

class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind),dimension(km3_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz)    &
                      ,intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
            end      subroutine boco_3d_gh  
              module          subroutine bocoT_2d_g1                          &
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
 import mg_intstate_type
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
              end      subroutine bocoT_2d_g1 
            module            subroutine bocoT_2d_gh                          &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------

 import mg_intstate_type
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
            end             subroutine bocoT_2d_gh

                 module       subroutine bocoTx_2d_g1                         &
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
implicit none

class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
                 end       subroutine bocoTx_2d_g1   
                 module       subroutine bocoTx_2d_gh                         &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
                 end       subroutine bocoTx_2d_gh       
!-----------------------------------------------------------------------

          module              subroutine bocoTy_2d_g1                         &
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
implicit none

class(mg_intstate_type),target::this
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
          end              subroutine bocoTy_2d_g1   

                      module  subroutine bocoTy_2d_gh                         &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
                    end    subroutine bocoTy_2d_gh 

             module           subroutine bocoT_3d_g1                          &
(this,W,km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in)
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3_in,im_in,jm_in,Lm_in,nbx,nby,nbz
real(r_kind), dimension(km3_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz)   &
                       ,intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
             end            subroutine bocoT_3d_g1     
             module           subroutine bocoT_3d_gh                          &
(this,W,km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,Fimax_in,Fjmax_in,mygen_min,mygen_max)

implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,Lm_in,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind), dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby,1-nbz:Lm_in+nbz)    &
                       ,intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
             end           subroutine bocoT_3d_gh
                module        subroutine bocox_2d_gh                          &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------

implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
                end        subroutine bocox_2d_gh     
                 module       subroutine bocox_2d_g1                          &
(this,W,km_in,im_in,jm_in,nbx,nby)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
                 end       subroutine bocox_2d_g1 

        module                subroutine bocoy_2d_g1                          &
(this,W,km_in,im_in,jm_in,nbx,nby)

implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
        end                 subroutine bocoy_2d_g1
                  module      subroutine bocoy_2d_gh                          &
(this,W,km_in,im_in,jm_in,nbx,nby,Fimax_in,Fjmax_in,mygen_min,mygen_max)
!-----------------------------------------------------------------------
 
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km_in,im_in,jm_in,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km_in,1-nbx:im_in+nbx,1-nby:jm_in+nby),intent(inout):: W
integer(i_kind), dimension(this%gm), intent(in):: Fimax_in,Fjmax_in
                  end      subroutine bocoy_2d_gh   


                   module     subroutine upsend_all_g1                        &
(this,Harray,Warray,km_in)
!-----------------------------------------------------------------------
 import mg_intstate_type
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(in):: Harray
real(r_kind), dimension(km_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: Warray
                   end     subroutine upsend_all_g1                      
             module           subroutine upsend_all_gh                        &
(this,Harray,Warray,km_in,mygen_dn,mygen_up)
 import mg_intstate_type
implicit none
class(mg_intstate_type),target::this

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(in):: Harray
real(r_kind), dimension(km_in,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: Warray
integer(i_kind),intent(in):: mygen_dn,mygen_up
            end         subroutine upsend_all_gh   

            module            subroutine downsend_all_gh                      &
(this,Warray,Harray,km_in,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm),intent(in):: Warray
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
            end   subroutine downsend_all_gh 
              module       subroutine downsend_all_g2                      &
!                                                                      *
(this,Warray,Harray,km_in)
!-----------------------------------------------------------------------
implicit none
class(mg_intstate_type),target::this
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm),intent(in):: Warray
real(r_kind), dimension(km_in,1:this%imL,1:this%jmL),intent(out):: Harray
              end        subroutine downsend_all_g2    

end interface
!clt from mg_filtering
interface    
  module subroutine mg_filtering_procedure (this,mg_filt)
   import mg_intstate_type
    class(mg_intstate_type),target::this
     integer(i_kind),intent(in):: mg_filt

  end subroutine mg_filtering_procedure 
  module subroutine mg_filtering_rad1(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this
  end subroutine mg_filtering_rad1 
  module subroutine mg_filtering_rad2(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_rad2
  module subroutine mg_filtering_rad3(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_rad3
  module subroutine mg_filtering_lin1(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_lin1
module subroutine mg_filtering_lin2(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_lin2
module subroutine mg_filtering_lin3(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_lin3
module   subroutine mg_filtering_fast(this)
   import mg_intstate_type
    class(mg_intstate_type),target::this

  end subroutine mg_filtering_fast
module   subroutine sup_vrbeta1  & 
 (this,kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
   import mg_intstate_type
implicit none
        class(mg_intstate_type),target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

  end subroutine sup_vrbeta1
module   subroutine sup_vrbeta1T                        &
(this,kmax,hx,hy,hz,im,jm,lm,  pasp,ss, V) 
   import mg_intstate_type
    class(mg_intstate_type),target::this
  integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss
    

  end subroutine sup_vrbeta1T
 module  subroutine sup_vrbeta3                        &
  (this,kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
   import mg_intstate_type
    class(mg_intstate_type),target::this
    integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,this%i0:im,this%j0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(this%i0:im,this%j0:jm,1:lm), intent(in):: ss


  end subroutine sup_vrbeta3 
  module    subroutine sup_vrbeta3T                       &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(this,kmax,hx,hy,hz,im,jm,lm, pasp,ss,V)
!----------------------------------------------------------------------
   import mg_intstate_type
implicit none
        class(mg_intstate_type),target::this

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,this%i0-hx:im+hx,this%j0-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,this%i0:im,this%j0:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(this%i0:im,this%j0:jm,1:lm), intent(in):: ss
 end subroutine sup_vrbeta3T 


 end interface
!clt from mg_generations.f90
  interface
    module                    subroutine upsending_all                        &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!                                                                      !
!***********************************************************************
(this,V,H,lquart)
!-----------------------------------------------------------------------
   import mg_intstate_type
class (mg_intstate_type),target:: this
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
logical, intent(in):: lquart
end subroutine upsending_all

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
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
logical, intent(in):: lquart
end subroutine downsending_all 
 module    subroutine weighting_all                        &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H,lhelm)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
logical, intent(in):: lhelm
end subroutine  weighting_all
 module  subroutine upsending                            &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H

real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: V_INT
real(r_kind),dimension(this%km,this%i0-2:this%imL+2,this%j0-2:this%jmL+2):: H_INT
end subroutine upsending
 module   subroutine downsending                          &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(this,H,V)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V

end subroutine downsending
 module subroutine upsending2                           &
!***********************************************************************
!                                                                      !
!  Adjoint interpolate and upsend:                                     !
!       First from g1->g2 (V -> H)                                     !
!       Then  from g2->...->gn  (H -> H)                               !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(in):: V
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(out):: H
end subroutine upsending2

 module subroutine downsending2                         &
!***********************************************************************
!                                                                      !
!  Downsend, interpolate and add:                                      !
!      First from gm->g3...->g2                                        !
!      Then  from g2->g1                                               !
!                                                                      !
!***********************************************************************
(this,H,V)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(inout):: H
real(r_kind),dimension(this%km,1-this%hx:this%im+this%hx,1-this%hy:this%jm+this%hy),intent(inout):: V
end subroutine downsending2
module  subroutine weighting_helm                       &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
end subroutine weighting_helm 



 module subroutine weighting                            &
!***********************************************************************
!                                                                      !
!  Apply 2D differential operator to compound variable                 !
!                                                                      !
!***********************************************************************
(this,V,H)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this

real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: V
real(r_kind),dimension(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy),intent(inout):: H
end subroutine weighting

module subroutine adjoint                              &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,F,W,km_in,g)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,this%i0:this%im,this%j0:this%jm), intent(in):: F
real(r_kind), dimension(km_in,this%i0-2:this%imL+2,this%j0-2:this%jmL+2), intent(out):: W
end subroutine adjoint

module subroutine direct1                              &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using linearly squared interpolations                              !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,W,F,km_in,g)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,this%i0-2:this%imL+2,this%j0-2:this%jmL+2), intent(in):: W
real(r_kind), dimension(km_in,this%i0:this%im,this%j0:this%jm), intent(out):: F
end subroutine direct1
module subroutine adjoint2                             &
!***********************************************************************
!                                                                      !
!   Mapping from the high to low resolution grid                       !
!   using quadratics interpolations                                    !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,F,W,km_in,g)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,1:this%im,1:this%jm), intent(in):: F
real(r_kind), dimension(km_in,0:this%imL+1,0:this%jmL+1), intent(out):: W
end subroutine adjoint2

      module                  subroutine direct2                              &
!***********************************************************************
!                                                                      !
!   Mapping from the low to high resolution grid                       !
!   using quadratic interpolations                                     !
!                         - offset version -                           !
!                                                                      !
!***********************************************************************
(this,W,F,km_in,g)
!-----------------------------------------------------------------------
   import mg_intstate_type
implicit none
class (mg_intstate_type),target:: this
integer(i_kind),intent(in):: g
integer(i_kind),intent(in):: km_in
real(r_kind), dimension(km_in,0:this%imL+1,0:this%jmL+1), intent(in):: W
real(r_kind), dimension(km_in,1:this%im,1:this%jm), intent(out):: F
end subroutine direct2
  
               module         subroutine mg_initialize(this,inputfilename,obj_parameter)
   import mg_intstate_type
   import mg_parameter_type
class (mg_intstate_type):: this
character*(*),optional,intent(in) :: inputfilename
class(mg_parameter_type),optional,intent(in)::obj_parameter
   end subroutine mg_initialize
  module      subroutine mg_finalize(this)
   import mg_intstate_type
implicit none
class (mg_intstate_type)::this
        end subroutine mg_finalize
   



  end interface 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!cltthink doublecheck                        subroutine allocate_mg_intstate(this,km)
                        subroutine allocate_mg_intstate(this)
!***********************************************************************
!                                                                      !
! Allocate internal state variables                                    !
!                                                                      !
!***********************************************************************
   import mg_intstate_type
implicit none
class(mg_intstate_type),target::this


allocate(this%V(this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy,this%lm)) ; this%V=0.
allocate(this%VALL(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))  ; this%VALL=0.
allocate(this%HALL(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))  ; this%HALL=0.


allocate(this%a_diff_f(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))  ; this%a_diff_f=0. 
allocate(this%a_diff_h(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))  ; this%a_diff_h=0. 
allocate(this%b_diff_f(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy)) ; this%b_diff_f=0. 
allocate(this%b_diff_h(this%km,this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy)) ; this%b_diff_h=0. 

allocate(this%p_eps(this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))    ; this%p_eps=0.
allocate(this%p_del(this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))    ; this%p_del=0.
allocate(this%p_sig(this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))    ; this%p_sig=0.
allocate(this%p_rho(this%i0-this%hx:this%im+this%hx,this%j0-this%hy:this%jm+this%hy))    ; this%p_rho=0.

allocate(this%paspx(1,1,this%i0:this%im))                       ; this%paspx=0.
allocate(this%paspy(1,1,this%j0:this%jm))                                          ; this%paspy=0.

allocate(this%pasp1(1,1,1:this%lm))                                           ; this%pasp1=0.
allocate(this%pasp2(2,2,this%i0:this%im,this%j0:this%jm))                                    ; this%pasp2=0.
allocate(this%pasp3(3,3,this%i0:this%im,this%j0:this%jm,1:this%lm))               ; this%pasp3=0.

allocate(this%vpasp2(0:2,this%i0:this%im,this%j0:this%jm))                                   ; this%vpasp2=0.
allocate(this%hss2(this%i0:this%im,this%j0:this%jm,1:3))                                     ; this%hss2= 0.

allocate(this%vpasp3(1:6,this%i0:this%im,this%j0:this%jm,1:this%lm))                              ; this%vpasp3= 0.
allocate(this%hss3(this%i0:this%im,this%j0:this%jm,1:this%lm,1:6))                                ; this%hss3= 0.

allocate(this%ssx(this%i0:this%im))                                             ; this%ssx=0.
allocate(this%ssy(this%j0:this%jm))                                             ; this%ssy=0.
allocate(this%ss1(1:this%lm))                                             ; this%ss1=0.
allocate(this%ss2(this%i0:this%im,this%j0:this%jm))                                        ; this%ss2=0.
allocate(this%ss3(this%i0:this%im,this%j0:this%jm,1:this%lm))                                   ; this%ss3=0.

allocate(this%dixs(this%i0:this%im,this%j0:this%jm,3))                                     ; this%dixs=0
allocate(this%diys(this%i0:this%im,this%j0:this%jm,3))                                     ; this%diys=0

allocate(this%dixs3(this%i0:this%im,this%j0:this%jm,1:this%lm,6))                               ; this%dixs3=0
allocate(this%diys3(this%i0:this%im,this%j0:this%jm,1:this%lm,6))                               ; this%diys3=0
allocate(this%dizs3(this%i0:this%im,this%j0:this%jm,1:this%lm,6))                               ; this%dizs3=0

allocate(this%qcols(0:7,this%i0:this%im,this%j0:this%jm,1:this%lm))                             ; this%qcols=0

!
! In stnadalone version
!
!allocate(r_vol(km,0:nm,0:mm,2))                             ; r_vol=0.
!
! ... but in global version there will be 
!     r_vol2 and r_vol3 for 2d and 3d variables
! and r_vol3 will need to be given vertical structure
!

!
!cltallocate(WORKA(km,n0:nm,m0:mm))                               ; WORKA=0.

!
! for re-decomposition
!

allocate(this%iref(this%n0:this%nm))                                     ; this%iref=0
allocate(this%jref(this%m0:this%mm))                                     ; this%jref=0

allocate(this%cx0(this%n0:this%nm))                                      ; this%cx0=0.
allocate(this%cx1(this%n0:this%nm))                                      ; this%cx1=0.
allocate(this%cx2(this%n0:this%nm))                                      ; this%cx2=0.
allocate(this%cx3(this%n0:this%nm))                                      ; this%cx3=0.

allocate(this%cy0(this%m0:this%mm))                                      ; this%cy0=0.
allocate(this%cy1(this%m0:this%mm))                                      ; this%cy1=0.
allocate(this%cy2(this%m0:this%mm))                                      ; this%cy2=0.
allocate(this%cy3(this%m0:this%mm))                                      ; this%cy3=0.

!TEST
!       call finishMPI
!TEST

allocate(this%p_coef(4))                                      ; this%p_coef=0.
allocate(this%q_coef(4))                                      ; this%q_coef=0.

allocate(this%a_coef(3))                                      ; this%a_coef=0.
allocate(this%b_coef(3))                                      ; this%b_coef=0.


allocate(this%cf00(this%n0:this%nm,this%m0:this%mm))                            ; this%cf00=0.
allocate(this%cf01(this%n0:this%nm,this%m0:this%mm))                            ; this%cf01=0.
allocate(this%cf02(this%n0:this%nm,this%m0:this%mm))                            ; this%cf02=0.
allocate(this%cf03(this%n0:this%nm,this%m0:this%mm))                            ; this%cf03=0.
allocate(this%cf10(this%n0:this%nm,this%m0:this%mm))                            ; this%cf10=0.
allocate(this%cf11(this%n0:this%nm,this%m0:this%mm))                            ; this%cf11=0.
allocate(this%cf12(this%n0:this%nm,this%m0:this%mm))                            ; this%cf12=0.
allocate(this%cf13(this%n0:this%nm,this%m0:this%mm))                            ; this%cf13=0.
allocate(this%cf20(this%n0:this%nm,this%m0:this%mm))                            ; this%cf20=0.
allocate(this%cf21(this%n0:this%nm,this%m0:this%mm))                            ; this%cf21=0.
allocate(this%cf22(this%n0:this%nm,this%m0:this%mm))                            ; this%cf22=0.
allocate(this%cf23(this%n0:this%nm,this%m0:this%mm))                            ; this%cf23=0.
allocate(this%cf30(this%n0:this%nm,this%m0:this%mm))                            ; this%cf30=0.
allocate(this%cf31(this%n0:this%nm,this%m0:this%mm))                            ; this%cf31=0.
allocate(this%cf32(this%n0:this%nm,this%m0:this%mm))                            ; this%cf32=0.
allocate(this%cf33(this%n0:this%nm,this%m0:this%mm))                            ; this%cf33=0.

allocate(this%Lref(1:this%lm))                                   ; this%Lref=0
allocate(this%Lref_h(1:this%lm))                                 ; this%Lref_h=0

allocate(this%cvf1(1:this%lm))                                   ; this%cvf1=0.
allocate(this%cvf2(1:this%lm))                                   ; this%cvf2=0.
allocate(this%cvf3(1:this%lm))                                   ; this%cvf3=0.
allocate(this%cvf4(1:this%lm))                                   ; this%cvf4=0.

allocate(this%cvh1(1:this%lm))                                   ; this%cvh1=0.
allocate(this%cvh2(1:this%lm))                                   ; this%cvh2=0.
allocate(this%cvh3(1:this%lm))                                   ; this%cvh3=0.
allocate(this%cvh4(1:this%lm))                                   ; this%cvh4=0.


!-----------------------------------------------------------------------
                        endsubroutine allocate_mg_intstate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_mg_weights(this)
!***********************************************************************
!                                                                      !
! Define weights and scales                                            !
   import mg_intstate_type
!                                                                      !
implicit none
class (mg_intstate_type),target::this
!***********************************************************************
integer(i_kind):: i,j,L
real(r_kind):: gen_fac
!-----------------------------------------------------------------------

      this%p_eps(:,:)=0.0
      this%p_del(:,:)=0.0
      this%p_sig(:,:)=0.0
      this%p_rho(:,:)=0.0

!--------------------------------------------------------
      gen_fac=1.
      this%a_diff_f(:,:,:)=this%mg_weig1 
      this%a_diff_h(:,:,:)=this%mg_weig1 

      this%b_diff_f(:,:,:)=0.
      this%b_diff_h(:,:,:)=0.

!      r_vol(:,:,:,1)=1.


      select case(this%my_hgen)
        case(2) 
!          r_vol(:,:,:,2)=0.25             ! In standalone case
!          gen_fac=0.25
          this%a_diff_h(:,:,:)=this%mg_weig2
          this%b_diff_h(:,:,:)=0.
        case(3) 
!          r_vol(:,:,:,2)=0.0625           ! In standalone case
!          gen_fac=0.0625
          this%a_diff_h(:,:,:)=this%mg_weig3 
          this%b_diff_h(:,:,:)=0.
        case default 
!          r_vol(:,:,:,2)=0.015625         ! In standalone case
!          gen_fac=0.015625
          this%a_diff_h(:,:,:)=this%mg_weig4
          this%b_diff_h(:,:,:)=0.
      end select


          do L=1,this%lm
           this%pasp1(1,1,L)=this%pasp01
          enddo

          do i=this%i0,this%im
            this%paspx(1,1,i)=this%pasp02
          enddo  
          do j=this%j0,this%jm
            this%paspy(1,1,j)=this%pasp02
          enddo  

          do j=this%i0,this%jm
          do i=this%j0,this%im
            this%pasp2(1,1,i,j)=this%pasp02*(1.+this%p_del(i,j))
            this%pasp2(2,2,i,j)=this%pasp02*(1.-this%p_del(i,j))
            this%pasp2(1,2,i,j)=this%pasp02*this%p_eps(i,j)     
            this%pasp2(2,1,i,j)=this%pasp02*this%p_eps(i,j)     
          end do
          end do

        do L=1,this%lm
          do j=this%i0,this%jm
          do i=this%j0,this%im
            this%pasp3(1,1,i,j,l)=this%pasp03*(1+this%p_del(i,j))
            this%pasp3(2,2,i,j,l)=this%pasp03
            this%pasp3(3,3,i,j,l)=this%pasp03*(1-this%p_del(i,j))
            this%pasp3(1,2,i,j,l)=this%pasp03*this%p_eps(i,j)
            this%pasp3(2,1,i,j,l)=this%pasp03*this%p_eps(i,j)
            this%pasp3(2,3,i,j,l)=this%pasp03*this%p_sig(i,j)
            this%pasp3(3,2,i,j,l)=this%pasp03*this%p_sig(i,j)
            this%pasp3(1,3,i,j,l)=this%pasp03*this%p_rho(i,j)
            this%pasp3(3,1,i,j,l)=this%pasp03*this%p_rho(i,j)
          end do
          end do
        end do


        call this%cholaspect(1,this%lm,this%pasp1)
        call this%cholaspect(this%i0,this%im,this%j0,this%jm,this%pasp2)
        call this%cholaspect(this%i0,this%im,this%j0,this%jm,1,this%lm,this%pasp3)


        call this%getlinesum(this%hx,this%i0,this%im,this%paspx,this%ssx)
        call this%getlinesum(this%hy,this%j0,this%jm,this%paspy,this%ssy)
        call this%getlinesum(this%hz,1,this%lm,this%pasp1,this%ss1)
        call this%getlinesum(this%hx,this%i0,this%im,this%hy,this%j0,this%jm,this%pasp2,this%ss2)
        call this%getlinesum(this%hx,this%i0,this%im,this%hy,this%j0,this%jm,this%hz,1,this%lm,this%pasp3,this%ss3)
!-----------------------------------------------------------------------
                        endsubroutine def_mg_weights

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_line(this)
   import mg_intstate_type
implicit none
class(mg_intstate_type),target::this
integer(i_kind):: i,j,L,icol
logical:: ff
!***********************************************************************
!                                                                      !
! Inititate line filters                                               !
!                                                                      !
!***********************************************************************
!-----------------------------------------------------------------------

  do j=this%j0,this%jm
  do i=this%i0,this%im
    call t22_to_3(this%pasp2(:,:,i,j),this%vpasp2(:,i,j))
  enddo
  enddo

  do l=1,this%lm
  do j=this%j0,this%jm
  do i=this%i0,this%im
    call t33_to_6(this%pasp3(:,:,i,j,l),this%vpasp3(:,i,j,l))
  enddo
  enddo
  enddo



  call inimomtab(this%p,this%nh,ff)

  call tritform(this%i0,this%im,this%i0,this%jm,this%vpasp2, this%dixs,this%diys, ff)

  do icol=1,3
    this%hss2(:,:,icol)=this%vpasp2(icol-1,:,:)
  enddo  


  call hextform(this%i0,this%im,this%j0,this%jm,1,this%lm,this%vpasp3,this%qcols,this%dixs3,this%diys3,this%dizs3, ff)


  do icol=1,6
    this%hss3(:,:,:,icol)=this%vpasp3(icol,:,:,:)
  enddo
 

!-----------------------------------------------------------------------
                        endsubroutine init_mg_line

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine deallocate_mg_intstate(this)
implicit none
class (mg_intstate_type),target:: this
!***********************************************************************
!                                                                      !
! Deallocate internal state variables                                  !
!                                                                      !
!***********************************************************************

deallocate(this%V)

deallocate(this%HALL,this%VALL)

deallocate(this%a_diff_f,this%b_diff_f)
deallocate(this%a_diff_h,this%b_diff_h)
deallocate(this%p_eps,this%p_del,this%p_sig,this%p_rho,this%pasp1,this%pasp2,this%pasp3,this%ss1,this%ss2,this%ss3)
deallocate(this%dixs,this%diys)
deallocate(this%dixs3,this%diys3,this%dizs3)
deallocate(this%qcols)
!
! for testing
!
!cltthink deallocate(WORKA)

!
! for re-decomposition
!
deallocate(this%iref,this%jref)

deallocate(this%cf00,this%cf01,this%cf02,this%cf03,this%cf10,this%cf11,this%cf12,this%cf13)
deallocate(this%cf20,this%cf21,this%cf22,this%cf23,this%cf30,this%cf31,this%cf32,this%cf33)

deallocate(this%Lref,this%Lref_h)

deallocate(this%cvf1,this%cvf2,this%cvf3,this%cvf4)

deallocate(this%cvh1,this%cvh2,this%cvh3,this%cvh4)

deallocate(this%cx0,this%cx1,this%cx2,this%cx3)
deallocate(this%cy0,this%cy1,this%cy2,this%cy3)

deallocate(this%p_coef,this%q_coef)
deallocate(this%a_coef,this%b_coef)



                        end subroutine deallocate_mg_intstate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        end module mg_intstate
