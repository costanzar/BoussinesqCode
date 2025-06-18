!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module generates Initial Conditions (IC) for vorticity and temperature (buoyancy) fields  !!
!! It can read a NetCDF as IC                                                                     !!
!! or generated IC based on a description in physical or Fourier space                            !!
!! or IC can be a superposition of the above                                                      !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module init_condition
!  use param
  use param_fftw
  use IO_netcdf
!  use velvorproj
  use nm_decomp
  implicit none
  ! -----------------------------------------------------------------  
  ! PARAMETERS (icmode set in param.F90)
  ! >>>> parameters for peak in energy spectrum
  ! ki is the index of wavenumber at which IC peaks (in energy spc) ONLY for icmode = 0
  real,    parameter :: ki     =   8.0
  real,    parameter :: ke_i   =   2.0     !  initial kinetic energy
  
  ! Make internal variables and functions private
  PRIVATE :: ranno, ran1

CONTAINS
  
subroutine init_cond(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,zxr,zyr,zzr,ttr,ur,vr,wr,ts)
! Define the initial condition
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,uu,vv,ww,ge,g1,g2
  real,    intent(out), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,ttr,ur,vr,wr
  real,    intent(out) :: ts
  complex, dimension(:,:,:), allocatable :: zx0, zy0, zz0, tt0
  

  if (icmode==2) then
     print*,'This set of scripts cannot handle two-part ICs'
     call MPI_Abort(mpi_comm_world,istatus)
!!$     ! superpostion of input NetCDF file and user defined fields
!!$     ! for user defined fields you can use the functions in this module like single_planewave ...
!!$     allocate(zx0(iktx,ikty,iktzp))
!!$     allocate(zy0(iktx,ikty,iktzp))
!!$     allocate(zz0(iktx,ikty,iktzp))
!!$     allocate(tt0(iktx,ikty,iktzp))
!!$     call ncreadin(zx,zy,zz,tt,ts)
!!$
!!$     ! superposing the wave field and background QG
!!$     zx = zx + zx0
!!$     zy = zy + zy0
!!$     zz = zz + zz0
!!$     tt = tt + tt0
!!$
!!$     ! free the memory for auxilliary variables 
!!$     deallocate(zx0,zy0,zz0,tt0)

  elseif (icmode==1) then
     call ncreadin(zx,zy,zz,tt,ts)
     if (mype.eq.0) then
        write(iuRESULT,*) '                '
        write(iuRESULT,*) 'Initial Condition: ---------------------------------------'
        write(iuRESULT,*) 'Starting from Zk.in.ncf'
        call flush(iuRESULT)
     endif
  else ! (icmode==0)
     call exp_peak(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2)
     
     if (mype.eq.0) then
        write(iuRESULT,*) '                '
        write(iuRESULT,*) 'Initial Condition: ---------------------------------------'
        write(iuRESULT,*) 'Starting from a Guassian Peak in Eng. Spec'
        write(iuRESULT,*) '            peak at ki = ', ki
        write(iuRESULT,*) '      initial kin eng. = ', ke_i
        call flush(iuRESULT)
     endif
  endif

  call wtoab(zx,zy,zz,tt,ge,g1,g2,uu,vv,ww)
  call proj(zx,zy,zz)
  
  return
end subroutine init_cond

subroutine exp_peak(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2)
!!! Initialize in Fourier space with Gausian peak at ki (defined above)  
  implicit none
  include 'mpif.h'
  complex, intent(out), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,uu,vv,ww,ge,g1,g2
  integer :: ikx,iky,ikz,ikza
  real    :: kx,ky,kz,wk,kh,khn,wkn,kzn
  real    :: phase,EK
  real    :: kinen,poten,vorten,waven
  
  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           wk  = sqrt(kx*kx + ky*ky + kz*kz)
           khn = sqrt((kx*L1/twopi)**2 + (ky*L2/twopi)**2)
           wkn = sqrt((kx*L1/twopi)**2 + (ky*L2/twopi)**2 + (kz*L3/twopi)**2)
           kzn = kz * L3/twopi
           if (L(ikx,iky,ikz).eq.1) then
              EK = exp(-((wkn-ki)/5.)**2)
              phase = ranno(0)
              zx(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              zy(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              zz(ikx,iky,ikz) = wk*sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              tt(ikx,iky,ikz) = sqrt(EK)*cexp(ZI*phase)
              phase = ranno(0)
              if (kx.eq.0.) then
                 zx(ikx,iky,ikz) = 1e-2*cmplx(real(zx(ikx,iky,ikz)),0.)
                 zy(ikx,iky,ikz) = 1e-2*cmplx(real(zy(ikx,iky,ikz)),0.)
                 zz(ikx,iky,ikz) = 1e-2*cmplx(real(zz(ikx,iky,ikz)),0.)
              endif
           endif
        enddo
     enddo
  enddo

  call wtoab(zx,zy,zz,tt,ge,g1,g2,uu,vv,ww)

  g1 = cmplx(0.,0.)
  g2 = cmplx(0.,0.)

  call atowb(ge,g1,g2,zx,zy,zz,tt,uu,vv,ww)
  call proj(zx,zy,zz)
 
  ! normalizing the initial condition
  call energy_full(zx,zy,zz,tt,uu,vv,ww,ge,g1,g2,kinen,poten,vorten,waven)
  zx = sqrt(ke_i/kinen)*zx
  zy = sqrt(ke_i/kinen)*zy
  zz = sqrt(ke_i/kinen)*zz 
  tt = sqrt(ke_i/kinen)*tt
     
end subroutine exp_peak



subroutine energy_full(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,ke,pe,eg,ea)
! Computes total energy (KE, PE, VE, WE) 
  implicit none 
  include 'mpif.h'

  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,ge,g1,g2
  real,    intent(out) :: ke,pe,eg,ea

  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,wk,wkh2,wkh2n,kzn
  real :: vh,vzx,vzy,vzz,vz
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,tmp

  call velo(zx,zy,zz,ux,uy,uz)
  call wtoab(zx,zy,zz,tt,ge,g1,g2,ux,uy,uz)

  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea  = 0.

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz*L3/twopi
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wkh2n = wkh2*(L1/twopi)**2
           wk = kx*kx + ky*ky + kz*kz

           if (L(ikx,iky,ikz).eq.1) then
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)
              
              vzx = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              ke = ke + vzx/wk
              vzy = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              ke = ke + vzy/wk
              vzz = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              ke = ke + vzz/wk
              vh = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              pe = pe  + vh
            
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )
            
              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv+vzx+vzy
                 zero_kh_geo = zero_kh_geo+VH*aj/bj
              endif
            
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo+VZX+VZY
                 zero_kz_grv = zero_kz_grv+VZZ+VH*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vz = real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))
                 eg = eg  + vz/wkh2
                 vz = real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
                 vz = real(g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))
                 ea = ea  + vz/wkh2
              endif

           endif
        enddo
     enddo
  enddo

  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj
  eg = eg + zero_kz_geo + zero_kh_geo
  ea = ea + zero_kz_grv + zero_kh_grv

  call mpi_allreduce(ke,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  ke=tmp
  call mpi_allreduce(pe,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  pe=tmp
  call mpi_allreduce(eg,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  eg=tmp
  call mpi_allreduce(ea,tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,istatus);  ea=tmp

end subroutine energy_full



function ranno (i)

! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  implicit none
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

   
end module init_condition

   
