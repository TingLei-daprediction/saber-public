
!*** 
! apply the solution from
!*** Namelist parameters
!***
mg_ampl01=>this%mg_ampl01
mg_ampl02=>this%mg_ampl02
mg_ampl03=>this%mg_ampl03
mg_weig1=>this%mg_weig1

mg_weig2=>this%mg_weig2
mg_weig3=>this%mg_weig3
mg_weig4=>this%mg_weig4
mgbf_proc=>this%mgbf_proc
mgbf_line=>this%mgbf_line
nxPE=>this%nxPE 
nyPE=>this%nyPE 
im_filt=>this%im_filt  
jm_filt=>this%jm_filt  
lquart=>this%lquart  
lhelm=>this%lhelm  

!*** 
!*** Number of generations
!***
gm=>this%gm              

!*** 
!*** Horizontal resolution 
!***

!
! Original number of data on GSI analysis grid
!
nA_max0=>this%nA_max0  
mA_max0=>this%mA_max0 

!
! Global number of data on Analysis grid
!
nm0 =>this%nm0         
 mm0=>this%mm0         

!
! Number of PEs on Analysis grid
!
nxm  =>this%nxm           
mym =>this%mym            

!
! Number of data on local Analysis grid
!
nm =>this%nm          
mm =>this%mm         

!
! Number of data on global Filter grid
!
im00=>this%im00  
jm00=>this%jm00  

!
! Number of data on local  Filter grid
!
im=>this%im  
jm =>this%jm     

!
! Initial index on local  Filter grid
!
i0=>this%i0  
j0=>this%j0      
!
! Initial index on local analysis grid
!
n0=>this%n0  
m0=>this%m0      

!
! Halo on local Filter grid 
!
ib=>this%ib  
jb=>this%jb           

!
! Halo on local Analysis grid 
!
 nb=>this%nb  
 mb=>this%mb       


hx=>this%hx  
hy=>this%hy  
hz=>this%hz  
p=>this%p  
nh=>this%nh
nfil=>this%nfil  
pasp01=>this%pasp01 
pasp02=>this%pasp02
pasp03=>this%pasp03
pee2=>this%pee2
rmom2_1=>this%rmom2_1
rmom2_2=>this%rmom2_2
rmom2_3=>this%rmom2_3
rmom2_4=>this%rmom2_4


maxpe_fgen=>this%maxpe_fgen  
ixm=>this%ixm  
jym=>this%jym
nxy=>this%nxy  
im0=>this%im0
jm0=>this%jm0
Fimax=>this%Fimax 
Fjmax=>this%Fjmax  
FimaxL=>this%FimaxL  
FjmaxL=>this%FjmaxL  

npes_filt=>this%npes_filt  

maxpe_filt=>this%maxpe_filt  

 imL=>this%imL
 jmL=>this%jmL  
lm=>this%lm              ! number of vertical layers
 lm05=>this%lm05            ! half of vertical levels
km2_f=>this%km2_f           ! number of 2d variables for filtering
km3_f=>this%km3_f           ! number of 3d variables for filtering 
km2_e=>this%km2_e           ! number of 2d variables for ensemble
km3_e=>this%km3_e           ! number of 3d variables for ensemble
l_filt=>this%l_filt               ! logical flag for filtering or enseble
 lmf =>this%lmf            ! number of vertical levels for filtering (generation one)
lmh =>this%lmh            ! number of vertical levels for filtering (high generations)



lengthx=>this%lengthx
lengthy=>this%lengthy

xa0=>this%xa0
ya0=>this%ya0
xf0=>this%xf0
yf0=>this%yf0
dxf=>this%dxf
dyf=>this%dyf
dxa=>this%dxa
dya=>this%dya

npadx=>this%npadx           ! x padding on analysis grid
mpady=>this%mpady           ! y padding on analysis grid

ipadx=>this%ipadx           ! x padding on filter decomposition
jpady=>this%jpady           ! y padding on filter deocmposition

!
! Just for standalone test
!
!clt from entrymode.f90
ldelta=>this%ldelta  
km=>this%km 
km2=>this%km2
km3=>this%km3
!clt from mg_mpstuff.f90
mype=>this%mype  
c_mype=>this%c_mype  
npes=>this%npes 
iTYPE=>this%iTYPE
rTYPE=>this%rTYPE 
dTYPE=>this%dTYPE 
mpi_comm_comp=>this%mpi_comm_comp
ierror=>this%ierror  
mpi_comm_work=>this%mpi_comm_work 
 group_world=>this%group_world
 group_work=>this%group_work

 mype_gr=>this%mype_gr
 npes_gr=>this%npes_gr

 my_hgen=>this%my_hgen
 mype_hgen=>this%mype_hgen

 l_hgen=>this%l_hgen
 nx=>this%nx
 my=>this%my
!clt moved from *_mg_domain.f90
!logical,dimension(2):: Flwest,Fleast,Flnorth,Flsouth
   Flwest=>this%Flwest
  Fleast=>this%Fleast
  Flnorth=>this%Flnorth
  Flsouth=>this%Flsouth
!integer(i_kind),dimension(2):: Fitarg_n,Fitarg_e,Fitarg_s,Fitarg_w                         
  Fitarg_n=>this%Fitarg_n
  Fitarg_e=>this%Fitarg_e
  Fitarg_s=>this%Fitarg_s
  Fitarg_w=>this%Fitarg_w
!integer(i_kind),dimension(2):: Fitarg_sw,Fitarg_se,Fitarg_ne,Fitarg_nw
  Fitarg_sw=>this%Fitarg_sw
  Fitarg_se=>this%Fitarg_se
  Fitarg_ne=>this%Fitarg_ne
  Fitarg_nw=>this%Fitarg_nw

!logical,dimension(2):: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne
  Flsendup_sw=>this%Flsendup_sw
  Flsendup_se=>this%Flsendup_se
  Flsendup_nw=>this%Flsendup_nw
  Flsendup_ne=>this%Flsendup_ne
  Fitarg_up=>this%Fitarg_up

  itargdn_sw=>this%itargdn_sw
  itargdn_se=>this%itargdn_se
  itargdn_ne=>this%itargdn_ne
  itargdn_nw=>this%itargdn_nw


  itarg_wA=>this%itarg_wA
  itarg_eA=>this%itarg_eA
  itarg_sA=>this%itarg_sA
  itarg_nA=>this%itarg_nA
  lwestA=>this%lwestA
  leastA=>this%leastA
  lsouthA=>this%lsouthA
  lnorthA=>this%lnorthA


  ix=>this%ix
  jy=>this%jy


  mype_filt=>this%mype_filt

