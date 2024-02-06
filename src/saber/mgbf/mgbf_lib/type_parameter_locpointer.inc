! https://stackoverflow.com/questions/24990491/fortran-namelist-associate-does-not-work
real(r_kind),pointer :: mg_ampl01,mg_ampl02,mg_ampl03
real(r_kind),pointer:: mg_weig1,mg_weig2,mg_weig3,mg_weig4
integer(i_kind),pointer:: mgbf_proc
logical,pointer:: mgbf_line
integer(i_kind),pointer:: nxPE,nyPE,im_filt,jm_filt
logical,pointer:: lquart,lhelm
integer(i_kind),pointer:: gm            
integer(i_kind),pointer:: nA_max0
integer(i_kind),pointer:: mA_max0
integer(i_kind),pointer:: nm0        
integer(i_kind),pointer:: mm0       
integer(i_kind),pointer:: nxm           
integer(i_kind),pointer:: mym           
integer(i_kind),pointer:: nm         
integer(i_kind),pointer:: mm        
integer(i_kind),pointer:: im00
integer(i_kind),pointer:: jm00
integer(i_kind),pointer:: im
integer(i_kind),pointer:: jm    
integer(i_kind),pointer:: i0
integer(i_kind),pointer:: j0    
integer(i_kind),pointer:: n0
integer(i_kind),pointer:: m0    
integer(i_kind),pointer:: ib
integer(i_kind),pointer:: jb         
integer(i_kind),pointer:: nb
integer(i_kind),pointer:: mb     
integer(i_kind),pointer:: hx,hy,hz
integer(i_kind),pointer:: p
integer(i_kind),pointer:: nh,nfil
real(r_kind),pointer:: pasp01,pasp02,pasp03
real(r_kind),pointer:: pee2,rmom2_1,rmom2_2,rmom2_3,rmom2_4
integer, pointer, dimension(:):: maxpe_fgen
integer, pointer, dimension(:):: ixm,jym,nxy
integer, pointer, dimension(:):: im0,jm0
integer, pointer, dimension(:):: Fimax,Fjmax
integer, pointer, dimension(:):: FimaxL,FjmaxL
integer(i_kind),pointer:: npes_filt
integer(i_kind),pointer:: maxpe_filt
integer(i_kind),pointer:: imL,jmL
integer(i_kind),pointer:: lm            ! number of vertical layers
integer(i_kind),pointer:: lm05          ! half of vertical levels
integer(i_kind),pointer:: km2_f         ! number of 2d variables for filtering
integer(i_kind),pointer:: km3_f         ! number of 3d variables for filtering 
integer(i_kind),pointer:: km2_e         ! number of 2d variables for ensemble
integer(i_kind),pointer:: km3_e         ! number of 3d variables for ensemble
logical,pointer :: l_filt               ! logical flag for filtering or enseble
!integer(i_kind):: km            ! number of vertically stacked all variables (km=km2+lm*km3)
integer(i_kind),pointer:: lmf           ! number of vertical levels for filtering (generation one)
integer(i_kind),pointer:: lmh           ! number of vertical levels for filtering (high generations)
real(r_kind),pointer:: lengthx,lengthy,xa0,ya0,xf0,yf0
real(r_kind),pointer:: dxf,dyf,dxa,dya
integer(i_kind),pointer:: npadx         ! x padding on analysis grid
integer(i_kind),pointer:: mpady         ! y padding on analysis grid
integer(i_kind),pointer:: ipadx         ! x padding on filter decomposition
integer(i_kind),pointer:: jpady         ! y padding on filter deocmposition
logical,pointer:: ldelta
!clt from mg_entrymod.f90
integer(i_kind),pointer:: km,km2,km3
!clt from mg_mppstuff.f90
integer(i_kind),pointer:: mype
character(len=5),pointer:: c_mype
integer(i_kind),pointer:: npes,iTYPE,rTYPE,dTYPE,mpi_comm_comp,ierror
integer(i_kind),pointer:: mpi_comm_work,group_world,group_work

integer(i_kind),pointer:: mype_gr,npes_gr

integer(i_kind),pointer:: my_hgen
integer(i_kind),pointer:: mype_hgen
logical ,pointer:: l_hgen
integer(i_kind),pointer:: nx,my
!clt moved from *_mg_domain.f90
logical,dimension(:),pointer:: Flwest,Fleast,Flnorth,Flsouth
integer(i_kind),dimension(:),pointer:: Fitarg_n,Fitarg_e,Fitarg_s,Fitarg_w                         
integer(i_kind),dimension(:),pointer:: Fitarg_sw,Fitarg_se,Fitarg_ne,Fitarg_nw

logical,dimension(:),pointer:: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne
integer(i_kind),dimension(:),pointer:: Fitarg_up 

integer(i_kind),pointer:: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw


integer(i_kind),pointer:: itarg_wA,itarg_eA,itarg_sA,itarg_nA
logical,pointer:: lwestA,leastA,lsouthA,lnorthA


integer(i_kind),pointer:: ix,jy

integer(i_kind),dimension(:),pointer:: mype_filt


