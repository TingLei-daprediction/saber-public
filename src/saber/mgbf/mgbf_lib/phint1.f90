!#
!                                *********************************************
!                                *             MODULE phint1                 *
!                                *  R. J. Purser, NOAA/NCEP/EMC       2025   *
!                                *          jim.purser@noaa.gov              *
!                                *                                           *
!                                *********************************************
! Use interpolations of phint.f90 to construct a grid uniform in units of
! the "scale" given by a gridded profile. Also, use the interpolation of
! the logarithm of a given profile, followed by application of the exponential
! function, to ensure that the interpolation of a positive gridded function
! from one grid to another remains both smooth and positive.
!
! COMPILE AFTER: { phint.f90 }
!
!============================================================================
module phint1
!============================================================================
implicit none
private
public:: make_ssgrid, logintgrid

interface make_ssgrid
   module procedure make_ssgrid
end interface make_ssgrid
interface logintgrid
   module procedure logintgrid
end interface logintgrid

contains
   
!============================================================================
subroutine make_ssgrid(nz,nf,ns,sigofz,sstop,dss,isofz,zofis)!  [make_ssgrid]
!============================================================================
! Use the vertical profile, sigofs, of idealized correlation scale
! on the unit-spaced model grid of nz spaces to derive the total integrated
! depth, sstop, of the vertical domain in these scale units. Then,
! by using careful interpolations of the log of sigofz on the grid refined
! in the vertical by a factor of nf, divide the vertical domain into
! a new grid whose spacing is uniform in these scale units and which
! possesses ns grid spaces. On this grid, the correslation scale is
! constant, and can be taken to be sstopons=sstop/ns.
! Also, output the array, isofz, defining the index-coordinate of
! the new grid that corresponds to each model grid level, and the
! model grid index coordinate, zofis, that corresponds to each level of
! out new scale-grid. All grids are assumed to go from index 0.
!============================================================================
use pkind, only: dp,spi
use pietc, only: u1,o2
use phint, only: wint3,whint
implicit none
integer(spi),            intent(in ):: nz,nf,ns
real(dp),dimension(0:nz),intent(in ):: sigofz
real(dp),                intent(out):: sstop,dss
real(dp),dimension(0:nz),intent(out):: isofz
real(dp),dimension(0:ns),intent(out):: zofis
!----------------------------------------------------------------------------
real(dp),dimension(0:nz)   :: zs,logsig
real(dp),dimension(0:nz*nf):: zsf,logsigf,ssf
real(dp),dimension(0:ns)   :: ss
real(dp),dimension(3)      :: w3
real(dp),dimension(4)      :: w4
real(dp)                   :: r,s,z,dzf
integer(spi)               :: iz,izf,izfm,izfp,is,nzf
!============================================================================  
! Interpolate the log of the sigofz distribution to a finer grid:
do iz=0,nz
   zs(iz)=iz
   logsig(iz)=log(sigofz(iz))
enddo
dzf=u1/nf
nzf=nz*nf
do izf=0,nzf
   zsf(izf)=izf*dzf
enddo

do izf=0,nzf
   z=zsf(izf)
   iz=min(nz-1,max(0,floor(z)))
   if(iz==0)then
      call wint3(zs(0:2),z,w3)
      logsigf(izf)=dot_product(w3,logsig(0:2))
   elseif(iz==nz-1)then
      call wint3(zs(nz-2:nz),z,w3)
      logsigf(izf)=dot_product(w3,logsig(nz-2:nz))
   else
      call whint(zs(iz-1:iz+2),z,w4)
      logsigf(izf)=dot_product(w4,logsig(iz-1:iz+2))
   endif
enddo

ssf(0)=0
do izf=1,nzf
   izfm=izf-1
   ssf(izf)=ssf(izfm)+exp(-(logsigf(izfm)+logsigf(izf))*o2)*dzf
enddo
sstop   =ssf(nzf)

! define the new grid of ns spaces that uniformly divides the
! range of ss:
dss=sstop/ns
isofz(0)=0
isofz(nz)=ns
do iz=1,nz-1
   izf=iz*nf
   isofz(iz)=ssf(izf)/dss
enddo
do is=0,ns
   ss(is)=is*dss
enddo
zofis(0)=0
zofis(ns)=nz
izfp=1
do is=1,ns-1
   s=ss(is)
   do
      if(ssf(izfp)>=s)exit
      izfp=izfp+1
   enddo
   izf=izfp-1
   r=(s-ssf(izf))/(ssf(izfp)-ssf(izf))
   zofis(is)=(izf+r)/nf
enddo
end subroutine make_ssgrid

!============================================================================
subroutine logintgrid(nz,ns,zofs,az, as)!                        [logintgrid]
!============================================================================
! From a grid [0:nz] of positive values, az, use logarithms
! to ensure that the smooth interpolation to a new grid [0:ns]
! of target values, as, all remain positive. The array zofs
! defines the index z-grid coordinates of each of the s-grid points.
!============================================================================
use pkind, only: dp,spi
use phint, only: wint3,whint
implicit none
integer(spi),            intent(in ):: nz,ns
real(dp),dimension(0:ns),intent(in ):: zofs
real(dp),dimension(0:nz),intent(in ):: az
real(dp),dimension(0:ns),intent(out):: as
!----------------------------------------------------------------------------
real(dp),dimension(0:nz):: zs,logaz
real(dp),dimension(3)   :: w3! 3-point interpolation weights (at ends)
real(dp),dimension(4)   :: w4! 4-point interpolation weights (interior)
real(dp)                :: logas,z
integer(spi)            :: is,iz
!============================================================================
do iz=0,nz
   zs(iz)=iz
   logaz(iz)=log(az(iz))
enddo
do is=0,ns
   z=zofs(is)
   iz=min(nz-1,max(0,floor(z)))
   if(iz==0)then
      call wint3(zs(0:2),z,w3)
      logas=dot_product(w3,logaz(0:2))
   elseif(iz==nz-1)then
      call wint3(zs(nz-2:nz),z,w3)
      logas=dot_product(w3,logaz(nz-2:nz))
   else
      call whint(zs(iz-1:iz+2),z,w4)
      logas=dot_product(w4,logaz(iz-1:iz+2))
   endif
   as(is)=exp(logas)
enddo
end subroutine logintgrid

end module phint1
!#
