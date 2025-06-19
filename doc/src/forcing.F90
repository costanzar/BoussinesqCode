!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module contains the forcing and all related functions                                     !!
!! It can force the waves or geostrophic or all modes                                             !!
!! See Waite 2017 Physics of Fluids for more details on forcing scheme                            !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module forcing
  use param
  use velvorproj
  
  implicit none
  ! -----------------------------------------------------------------  
  ! PARAMETERS (for forcing)
  real,    parameter :: kf     =   3.                    ! wavenumber for forcing
  real,    parameter :: tau    =   900.                  ! time scale for forcing
  real,    parameter :: alpha  =   exp(-delt/tau)        ! memory parameter
  real,    parameter :: eps    =   4.e-8                 ! approx target dissipation rate for forcing amplitude
  real,    parameter :: ampv   =   1.2e-3/sqrt(tau)      ! random forcing amplitude
  real,    parameter :: ampw   =   0.                    ! random wave forcing amplitude (when used)
  integer, parameter :: nfmax  =   1200                  ! max number of forced modes
  
  ! Make internal variables and functions private
  ! PRIVATE :: 

CONTAINS

subroutine force(nzxk,nzyk,nzzk,nttk,gtau,nt)

! Use this after call convol.  Forcing is specified to be divergence-free
! Forcing follows Herring & Metais, 1989, J. Fluid Mech 202, pp. 97-115

  implicit none
  include 'mpif.h'
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: nzxk,nzyk,nzzk,nttk
  complex, intent(inout) :: gtau(1200,4)
  integer, intent(in)    :: nt
  
  integer :: ikx,iky,ikz,ikza,ik
  complex :: fu,fv,fw,ft,f1k,f2k,f3k,ftk
  real :: rang,beta,c,g,tf,theta
  real :: kx,ky,kz,wkh,wk,wkh2,sk
  real :: pk,pp,tmp

  beta  = sqrt(1.-alpha**2)
  c     = sqrt(bj/(2.*aj))
  ik    = 0
  pk    = 0.
  pp    = 0.

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           if (L(ikx,iky,ikz).eq.1) then 
              kx = kxa(ikx)
              wkh2 = kx**2 + ky**2
              wkh  = sqrt(wkh2) 
              wk   = sqrt(kx**2 + ky**2 + kz**2)
              theta = wkh/wk
              tf    = 1./sqrt2
       

              ! FORCE modes inside spherical shell centered on k=kf, AND exclude VSHF   
              if (abs(wkh*L1/twopi-KF).le.1. .and. abs(abs(kz)*L3/twopi-1.).lt.0.5) then
                 ik=ik+1                   
                 gtau(ik,1) = alpha*gtau(ik,1) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,2) = alpha*gtau(ik,2) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,3) = alpha*gtau(ik,3) + beta*cmplx(rang(0),rang(0))
                 gtau(ik,4) = alpha*gtau(ik,4) + beta*cmplx(rang(0),rang(0))
             
              !!! specify amplitudes of vort and wave forcing separately
                 G = 1.
                 Fv = ampv * G * gtau(ik,1)
                 Fw = ampw * G * ( gtau(ik,2) + gtau(ik,3))
                 Ft = ampw * G * (-gtau(ik,2) + gtau(ik,3))
                 sk = sqrt(BF2*wkh2 + cor2*kz**2)
                 
                 f1k  = - kx*kz*BF/sk*Fv - wk*ky/wkh/sqrt2*Fw + ZI*cor*kx*kz**2/sqrt2/sk/wkh*Ft
                 f2k  = - ky*kz*BF/sk*Fv + wk*kx/wkh/sqrt2*Fw + ZI*cor*ky*kz**2/sqrt2/sk/wkh*Ft
                 f3k  =   BF*wkh**2/sk*Fv - ZI*cor*wkh*kz/sqrt2/sk*Ft
                 ftk  = - C*sqrt2*ZI*cor*kz/sk*Fv + C*BF*wkh/sk*Ft
                 nzxk(ikx,iky,ikz) = nzxk(ikx,iky,ikz) + f1k
                 nzyk(ikx,iky,ikz) = nzyk(ikx,iky,ikz) + f2k
                 nzzk(ikx,iky,ikz) = nzzk(ikx,iky,ikz) + f3k
                 nttk(ikx,iky,ikz) = nttk(ikx,iky,ikz) + ftk

                  ! compute power from forcing
!!$                 pk = pk + real(conjg(zx(ikx,iky,ikz))*f1k/wk2)
!!$                 pk = pk + real(conjg(zy(ikx,iky,ikz))*f2k/wk2)
!!$                 pk = pk + real(conjg(zz(ikx,iky,ikz))*f3k/wk2)
!!$                 pp = pp + real(conjg(tt(ikx,iky,ikz))*ftk)*aj/bj


                 !!! new: force everything isotropically. Use ampv as global amplitude.
                 ! G = (wk*L1/twopi-(kf-1.))*(kf+1.-wk*L1/twopi) 
                 ! Fu = ampv * G * gtau(ik,1) 
                 ! Fv = ampv * G * gtau(ik,2) 
                 ! Fw = ampv * G * gtau(ik,3) 
                 ! Ft = ampv * G * gtau(ik,4) *sqrt(bj/aj)
                 ! nzxk(ikx,iky,ikz) = nzxk(ikx,iky,ikz) + zi*(ky*Fw-kz*Fv)
                 ! nzyk(ikx,iky,ikz) = nzyk(ikx,iky,ikz) + zi*(kz*Fu-kx*Fw)
                 ! nzzk(ikx,iky,ikz) = nzzk(ikx,iky,ikz) + zi*(kx*Fv-ky*Fu)
                 ! nttk(ikx,iky,ikz) = nttk(ikx,iky,ikz) + Ft
 

                 if (ik.eq.nfmax) then
                    print*,'Forcing too many modes!'
                    stop
                 endif
                 
              endif
           endif
        enddo
     enddo
  enddo
 
  pk=2.*pk
  pp=2.*pp

  ! for output times, write power to file
!!$  if (mod(nt,ndiagevery).eq. 0 ) then
!!$     call mpi_reduce(pk,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); pk=tmp
!!$     call mpi_reduce(pp,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus); pp=tmp
!!$     if (mype.eq.0) then
!!$        write(98,4999) time,pk,pp
!!$        call flush(98)
!!$     endif
!!$  endif
  return
  
4999 format(1x,f12.2,2x,2(e21.14,2x))
end subroutine force


function rang(i)

  ! If i ne 0 then initializes RANNO with seed = i
  ! If i eq 0 then draws a random GAUSSIAN number with 
  ! mean and std = 1

  implicit none
  include 'mpif.h'
  real :: rang
  integer :: i 
  real :: v1,v2,R,FAC,twopi,ranno
  external :: ranno

  twopi = 4.*asin(1.)
  
  if (i.ne.0) then
    v1 = ranno(i)
  else
200 v1 = 2.*(ranno(0)+twopi/2.)/twopi -1.
    v2 = 2.*(ranno(0)+twopi/2.)/twopi -1.
    r = v1**2. + v2**2.
    if (r.gt.1.) goto 200
    fac = sqrt( -2.*log(r)/r)
    rang = v1*fac
  endif
  return
end function rang


function ranno (i)

! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  implicit none
  include 'mpif.h'
  real :: ranno
  integer :: i
  integer :: junk,ihold
  real :: twopi,ran1
  save junk

  twopi = 4.*asin(1.)
  if (i.ne.0) then
    if (i.gt.0) i = - i
    junk  = i
    ranno = (ran1(i)-0.5)*twopi
  else
    junk  = junk - 1
    ihold = junk
    ranno = (ran1(ihold)-0.5)*twopi
  endif
  return
end function ranno


function ran1(idum)
  
  implicit none
  include 'mpif.h'
  real :: ran1
  integer :: idum

  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k,iv(32),iy
  save iv,iy
  data iv /ntab*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
    idum=max(-idum,1)
    do j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      if (j.le.ntab) iv(j)=idum
    enddo
    iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum.lt.0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
end function ran1


end module forcing

   
