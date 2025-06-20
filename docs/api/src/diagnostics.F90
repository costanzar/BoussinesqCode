!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!! This module contains the diagnostics such as eng time series, energy and transfer spectra    !!
!! as well as outputing run.list                                                                !!
!! One of the essential parts of the code that cannot be eliminated even in the lightest version!!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

module diagnostics
!  use param
!  use velvorproj
  use nm_decomp
  implicit none
  !! Wether or not keep eng.dat (time series of energys) or eps.dat (time series of dissipations)
  integer, parameter :: keepENG=1, keepEPS=0
  !! ID for energy time series (eng.dat) and dissipation (eps.dat)
  integer, parameter :: iuENG = 41, iuEPS = 42
  !! which ENERGY SPECTRA to keep:
  !! SPCH: horizontal energy spectra, SPCZ: vertical eng spectra, SPC: spherical eng spectra
  integer, parameter :: keepSPCH=1,keepSPCZ=1,keepSPC=0
  !! which TRANSFER SPECTRA to keep:
  integer, parameter :: keepTRNH=1,keepTRNZ=1,keepTRN=0
  !! A flag for decomposing the wave energy and transfers to potential and kinetic parts
  !! it is costly (both memory & computation) --> avoid it if not necessary for analysis
  integer, parameter :: potkinwv_flag = 0
  !! Triad Transfers: if triad transfers (GG->G, GA->G, AA-> G and ... ) are need turn it on!
  integer, parameter :: triads_on = 0
  
  !! the file IDs of each type of spectra
  integer, parameter :: iuTRNH=71, iuTRNZ=72, iuTRN=73
  !! the file IDs of each type of spectra
  integer, parameter :: iuSPCH=51, iuSPCZ=52, iuSPC=53
 
  
  ! Make internal variables and functions private
  PRIVATE :: keepSPCH,keepSPCZ,keepSPC,keepTRNH,keepTRNZ,keepTRN,potkinwv_flag,triads_on 
  !PRIVATE :: iuSPCH,iuSPCZ,iuSPC,iuTRNH,iuTRNZ,iuTRN

CONTAINS

subroutine prep_diagnostics()
  ! Open the binary (formatted) files for energy spectra, transfers, energy time series and ...
  implicit none
  include 'mpif.h'

  if (mype.eq.0) then
     if (keepENG ==1) open (iuENG,file='eng.dat',    form='formatted')
     if (keepEPS ==1) open (iuEPS,file='eps.dat',    form='formatted')
     if (keepSPCH==1) open (iuSPCH,file='spch.dat',  form='formatted')
     if (keepSPCZ==1) open (iuSPCZ,file='spcz.dat',  form='formatted')
     if (keepSPC ==1) open (iuSPC ,file='spc.dat' ,  form='formatted')
     if (keepTRNH==1) open (iuTRNH,file='trnh.dat',  form='formatted')
     if (keepTRNZ==1) open (iuTRNZ,file='trnz.dat',  form='formatted')
     if (keepTRN ==1) open (iuTRN ,file='trn.dat' ,  form='formatted')
  endif
  
  return
end subroutine prep_diagnostics

subroutine out_eng(zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,zxwv,zywv,zzwv,ttwv,nt)
  !! prints the energies and flow parameters in an output file (e.g. "run.list")
  !! also dumps the energies and dissipations as function of time in 'eng.dat' and 'eps.dat'
  !! Note: zx,zy,zz,tt,ge,g1,g2 would not change in this subroutine but
  !! ux,uy,uz,zxwv,zywv,zzwv,ttwv do change
  !! usually rhzx,rhzy,rhzz are used for zxwv,zywv,zzwv,ttwv to save memory
  !! Hence, pay attention when using "uk,vk,wk,nttk,rhzx,rhzy,rhzz"
  implicit none
  include 'mpif.h'

  integer, intent(in)  :: nt
  complex, intent(in), dimension(iktx,ikty,iktzp)    :: zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,zxwv,zywv,zzwv,ttwv

  complex, dimension(:,:,:), allocatable :: gewv
  integer :: ikx,iky,ikz,ikza
  real :: kx,ky,kz,wk2,wkh2,wkh2n,kzn
  real :: vh,vzx,vzy,vzz,zzx,zzy,zzz
  real :: rms_verv,rms_horv,rossby,fr_z,fr_h,ke,pe,e,eg,ea
  real :: epsk,epsp,eps,epskh,epskv,epsph,epspv,tmp
  real :: zero_kz_geo,zero_kz_grv,zero_kh_grv,zero_kh_geo,pewv,kewv,scaletime
  

  if (potkinwv_flag == 1) then
     allocate(gewv(iktx,ikty,iktzp))
  endif
  
  call velo(zx,zy,zz,ux,uy,uz)

  rms_verv = 0.
  rms_horv = 0.
  zero_kz_geo = 0.
  zero_kz_grv = 0.
  zero_kh_geo = 0.
  zero_kh_grv = 0.
  ke = 0.
  pe = 0.
  eg = 0.
  ea = 0.
  epsk = 0.
  epsp = 0.
  epskh = 0.
  epsph = 0.
  epskv = 0.
  epspv = 0.
  eps  = 0.
  
  do ikz = 1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kzn = kz * (L3/twopi)
     do iky = 1,ikty
        ky = kya(iky)
        do ikx = 1,iktx
           kx = kxa(ikx)
           wkh2 = kx*kx + ky*ky
           wk2 = kx*kx + ky*ky + kz*kz
           if (L(ikx,iky,ikz).eq.1) then
              wkh2n = wkh2 * (L1/twopi)**2

              zzx   = real(zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)))
              zzy   = real(zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)))
              zzz   = real(zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)))
              vzx   = real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)))
              vzy   = real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)))
              vzz   = real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)))
              vh    = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)))
              
              ke    = ke + (zzx+zzy+zzz)/wk2
              pe    = pe  + vh
              rms_verv=rms_verv+zzz
              rms_horv=rms_horv+zzx+zzy

              if (keepEPS ==1) then
                 epsk  = epsk + (visch*wkh2**ilap+viscz*kz**(2*ilap))*(vzx+vzy+vzz)
                 epskh = epskh + (visch*wkh2**ilap)*(vzx+vzy+vzz)
                 epskv = epskv + (viscz*kz**(2*ilap))*(vzx+vzy+vzz)
                 epsp  = epsp + (visch*wkh2**ilap+viscz*kz**(2*ilap))*vh*aj/bj
                 epsph = epsph + (visch*wkh2**ilap)*vh*aj/bj
                 epspv = epspv + (viscz*kz**(2*ilap))*vh*aj/bj
              endif
                           
              if(wkh2n.lt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 zero_kh_grv = zero_kh_grv + vzx + vzy
                 zero_kh_geo = zero_kh_geo + vh*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).lt.1.e-10) then
                 zero_kz_geo = zero_kz_geo + vzx + vzy
                 zero_kz_grv = zero_kz_grv + vzz + vh*aj/bj
              endif
              
              if(wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 eg=eg+real(ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)))/wkh2
                 ea=ea+real(g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz))+g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)))/wkh2
              endif              
           endif ! L
        enddo
     enddo
  enddo
  
  !! decomposing wave energy to a kinetic and potential part
  if (potkinwv_flag == 1) then 
     gewv=cmplx(0.,0.)
     call atowb(gewv,g1,g2,zxwv,zywv,zzwv,ttwv,ux,uy,uz)
     kewv = 0.
     pewv = 0. 
     do ikz = 1,iktzp
        ikza = mype*iktzp+ikz
        kz = kza(ikza)
        do iky = 1,ikty
           ky = kya(iky)
           do ikx = 1,iktx
              kx = kxa(ikx)
              wk2 = kx*kx + ky*ky + kz*kz
              if (L(ikx,iky,ikz).eq.1) then
                 kewv    = kewv + real( zxwv(ikx,iky,ikz)*conjg(zxwv(ikx,iky,ikz)))/wk2
                 kewv    = kewv + real( zywv(ikx,iky,ikz)*conjg(zywv(ikx,iky,ikz)))/wk2
                 kewv    = kewv + real( zzwv(ikx,iky,ikz)*conjg(zzwv(ikx,iky,ikz)))/wk2
                 pewv    = pewv  + real( ttwv(ikx,iky,ikz)*conjg(ttwv(ikx,iky,ikz)))
              endif
           enddo
        enddo
     enddo
  endif


  
  epsk   = 2.*epsk
  epsp   = 2.*epsp
  epskh  = 2.*epskh
  epsph  = 2.*epsph
  epskv  = 2.*epskv
  epspv  = 2.*epspv
  eg     = eg + zero_kz_geo + zero_kh_geo
  ea     = ea + zero_kz_grv + zero_kh_grv
  if (aj.ne.0. .and. bj.ne.0.) pe = aj*pe/bj
  if (aj.ne.0. .and. bj.ne.0.) pewv = aj*pewv/bj

  call mpi_reduce(ke,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          ke=tmp
  call mpi_reduce(pe,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          pe=tmp
  call mpi_reduce(eg,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          eg=tmp
  call mpi_reduce(ea,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);          ea=tmp
  call mpi_reduce(rms_verv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);    rms_verv=tmp
  call mpi_reduce(rms_horv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);    rms_horv=tmp
  if (keepEPS == 1) then
     call mpi_reduce(epsk,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        epsk=tmp
     call mpi_reduce(epsp,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        epsp=tmp
     call mpi_reduce(epskh,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epskh=tmp
     call mpi_reduce(epsph,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epsph=tmp
     call mpi_reduce(epskv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epskv=tmp
     call mpi_reduce(epspv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);       epspv=tmp
  endif
  if (potkinwv_flag == 1) then 
     deallocate(gewv)
     call mpi_reduce(kewv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        kewv=tmp
     call mpi_reduce(pewv,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus);        pewv=tmp
  endif

  if (mype.eq.0) then ! prep for output

    if (cor.gt.1.e-10) then
       rossby = sqrt(2.*rms_verv/cor2)
    else
       rossby = - 999.
    endif
    if (bf.gt.1.e-8) then
      fr_z = sqrt(rms_horv)/bf
      fr_h = sqrt(2.*rms_verv)/bf
    else
      fr_z = - 999.
      fr_h = - 999.
    endif
    eps    = epsk + epsp
    if (aj.ne.0. .and. bj.ne.0.) then
      e  = (pe + ke)
    else
      e = - 999.
   endif
   
   if (nt.eq.0) then
      if (potkinwv_flag == 1) then
         write(iuRESULT,5043)
         write(iuRESULT,5042) repeat('----------',15)
      else
         write(iuRESULT,5049)
         write(iuRESULT,5042) repeat('----------',13)
      endif
   endif

   scaletime = 3600.0

   if (potkinwv_flag == 1) then
      write(iuRESULT,5044) time/scaletime,ke,pe,kewv,pewv,e,eg,ea,rossby,fr_z,fr_h
      if (keepENG ==1) write(iuENG,5045) time/scaletime,ke,pe,kewv,pewv,e,eg,ea,rossby,fr_z,fr_h
   else
      write(iuRESULT,5048) time/scaletime,ke,pe,e,eg,ea,rossby,fr_z,fr_h
      if (keepENG ==1) write(iuENG,5046) time/scaletime,ke,pe,e,eg,ea,rossby,fr_z,fr_h
   endif
   call flush(iuRESULT)
   if (keepEPS ==1) write(iuEPS,5046) time/scaletime, epsk, epsp, eps, epskh, epskv, epsph, epspv
   if (keepENG ==1) call flush(iuENG)
   if (keepEPS ==1) call flush(iuEPS)
endif
  return

5042  format(1x,a91)
5043  format(8x,' T',8x,'KE',8x,'PE',8x,'KEWV',6x,'PEWV',6x,'E',9x,'GE',8x,'AE',8x,'Ro',8x,'Fr_z',7x,'Fr_h')
5049  format(8x,' T',8x,'KE',8x,'PE',8x,'E',9x,'GE',8x,'AE',8x,'Ro',8x,'Fr_z',7x,'Fr_h')  
5044  format(1x,f12.2,2x,10(f8.3,2x))
5048  format(1x,f12.2,2x,8(f8.3,2x))  
5045  format(1x,f12.2,2x,10(e11.4,1x))
5046  format(1x,f12.2,2x,8(e11.4,1x))
5047  format(1x,f12.2,2x,7(e21.14,1x))

end subroutine out_eng


subroutine spec (zx,zy,zz,tt,ux,uy,uz,ge,g1,g2,zxwv,zywv,zzwv,ttwv,ispec,iu)
! Calculates spectra using: 
! ispec = 1: total wavenumber k
! ispec = 2: horizontal wavenumber kh
! ispec = 3: vertical wavenumber kz
  
  implicit none 
  include 'mpif.h'

  integer, intent(in) :: ispec,iu
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,ge,g1,g2
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ux,uy,uz,zxwv,zywv,zzwv,ttwv
  
  integer :: ikx,iky,ikz,ikza,j,i,j0
  integer :: nspz,nspz13,nspz15
  real    :: kx,ky,kz,wk,wkh2,wkh,wkh2n,kzn,kz2
  real    :: vt,kvisc,vzx,vzy,vzz
  complex :: div
  real, dimension(:,:), allocatable      :: spz,spztot
  integer, dimension(0:kts)              :: n,ntot
  complex, dimension(:,:,:), allocatable ::  gewv

  if (potkinwv_flag == 1) then
     allocate(gewv(iktx,ikty,iktzp))
     allocate(spz(0:kts,15),spztot(0:kts,15))  
  else
     allocate(spz(0:kts,13),spztot(0:kts,13))
  endif
  
  nspz     = kts+1
  nspz13   = (kts+1)*13
  nspz15   = (kts+1)*15
  
  if (ispec.eq.1) then
    j0 = 1
  elseif (ispec.eq.2) then
    j0 = 0
  elseif (ispec.eq.3) then
    j0 = 0
  else
    print*,"ispec error"
    stop
  endif

  N     = 0
  spz   = 0.

  call velo (zx,zy,zz,ux,uy,uz)

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     kz2 = kz*kz
     kzn = kz*L3/twopi
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           wkh2  = kx*kx+ky*ky
           wkh2n = wkh2 * (L1/twopi)**2
           wkh   = sqrt(wkh2)
           wk    = sqrt(kx*kx+ky*ky+kz*kz)
           kvisc = visch*wkh2**ilap+viscz*kz2**ilap
           
           if (ispec.eq.1) then
              j = int(wk*L1/twopi+0.5)
           elseif (ispec.eq.2) then
              j = int(wkh*L1/twopi+0.5)
           elseif (ispec.eq.3) then
              j = int(abs(kz)*L3/twopi+0.5)
           endif
           
           if (L(ikx,iky,ikz).eq.1) then
              if (j.lt.j0 .or. j.gt.kts) print*,'SPEC: SCREW-UP.',j,'ktx= ',ktx,'ispec=',ispec
              wk   = max(wk,  1.e-15)
              wkh2 = max(wkh2,1.e-15)
              
              ! Kinetic and potential energy.
              vzx      = real( zx(ikx,iky,ikz)*conjg(zx(ikx,iky,ikz)) )
              vzy      = real( zy(ikx,iky,ikz)*conjg(zy(ikx,iky,ikz)) )
              vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )
              vt       = real( tt(ikx,iky,ikz)*conjg(tt(ikx,iky,ikz)) )
              
              spz(j,1) = spz(j,1) + vzx/wk**2 + vzy/wk**2 + vzz/wk**2
              spz(j,2) = spz(j,2) + vt*aj/bj
              
              ! KE and PE dissipation
              spz(j,7) = spz(j,7) + kvisc*(vzx/wk**2 + vzy/wk**2 + vzz/wk**2)     
              spz(j,8) = spz(j,8) + kvisc*vt*aj/bj
              

              ! Geo, ageo decompostition.
              ! k \in R_k
              if (wkh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
                 vzx      = real( ge(ikx,iky,ikz)*conjg(ge(ikx,iky,ikz)) )
                 vzy      = real( g1(ikx,iky,ikz)*conjg(g1(ikx,iky,ikz)) )
                 vzz      = real( g2(ikx,iky,ikz)*conjg(g2(ikx,iky,ikz)) )
                 
                 spz(j,4) = spz(j,4) + vzx/wkh2
                 spz(j,5) = spz(j,5) + vzy/wkh2 + vzz/wkh2
                 ! spz(j,7) = spz(j,7) + kvisc*vzx/wkh2
                 ! spz(j,8) = spz(j,8) + kvisc*vzy/wkh2 + kvisc*vzz/wkh2
              endif

              ! Special cases: i) k_h=0, ii) k_z=0.
              vzx=real(ux(ikx,iky,ikz)*conjg(ux(ikx,iky,ikz)) )
              vzy=real(uy(ikx,iky,ikz)*conjg(uy(ikx,iky,ikz)) )
              vzz=real(uz(ikx,iky,ikz)*conjg(uz(ikx,iky,ikz)) )

              ! k \in V_k
              if(wkh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
                 spz(j,4) = spz(j,4) + vzz + vt*aj/bj
                 spz(j,5) = spz(j,5) + vzx + vzy
                 ! spz(j,7) = spz(j,7) + kvisc*(vzz + vt*aj/bj)
                 ! spz(j,8) = spz(j,8) + kvisc*(vzx + vzy)     
              endif

              ! k \in B_k
              if(abs(kzn).lt.1.e-10.and.wkh2n.gt.1.e-10) then
                 spz(j,4) = spz(j,4) + vzx + vzy
                 spz(j,5) = spz(j,5) + vzz + vt*aj/bj
                 ! spz(j,7) = spz(j,7) + kvisc*(vzx + vzy)     
                 ! spz(j,8) = spz(j,8) + kvisc*(vzz + vt*aj/bj)
              endif
              
              ! Buoyancy Flux
              if (aj.gt.0.) then 
                 vt = aj*real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              else
                 vt = real(conjg(uz(ikx,iky,ikz))*tt(ikx,iky,ikz))
              endif
              spz(j,6) = spz(j,6) + vt
              
              ! KE decomposed into u/v/w
              spz(j,9)  = spz(j,9)  + vzx
              spz(j,10) = spz(j,10) + vzy
              spz(j,11) = spz(j,11) + vzz

              ! Rotational and divergent KE
              div      = kx*ux(ikx,iky,ikz)+ky*uy(ikx,iky,ikz)
              vzx      = real( div*conjg(div) )/wkh2
              vzz      = real( zz(ikx,iky,ikz)*conjg(zz(ikx,iky,ikz)) )/wkh2
              spz(j,12) = spz(j,12) + vzz
              spz(j,13) = spz(j,13) + vzx           
              
              n(j) = n(j) + 2 
              
           endif
        enddo
     enddo
  enddo
  
  spz(:,3) = spz(:,1) + spz(:,2)

  if (potkinwv_flag == 1) then 
     gewv=cmplx(0.,0.)
     call atowb(gewv,g1,g2,zxwv,zywv,zzwv,ttwv,ux,uy,uz)
     do ikz = 1,iktzp
        ikza = mype*iktzp+ikz
        kz = kza(ikza)
        do iky = 1,ikty
           ky = kya(iky)
           do ikx = 1,iktx
              kx = kxa(ikx)
              wk    = sqrt(kx*kx+ky*ky+kz*kz)
              wkh   = sqrt(kx*kx+ky*ky)
              if (ispec.eq.1) then
                 j = int(wk*L1/twopi+0.5)
              elseif (ispec.eq.2) then
                 j = int(wkh*L1/twopi+0.5)
              elseif (ispec.eq.3) then
                 j = int(abs(kz)*L3/twopi+0.5)
              endif
              
              if (L(ikx,iky,ikz).eq.1) then
                 ! Kinetic and potential energy of the waves
                 vzx      = real( zxwv(ikx,iky,ikz)*conjg(zxwv(ikx,iky,ikz)) )
                 vzy      = real( zywv(ikx,iky,ikz)*conjg(zywv(ikx,iky,ikz)) )
                 vzz      = real( zzwv(ikx,iky,ikz)*conjg(zzwv(ikx,iky,ikz)) )
                 vt       = real( ttwv(ikx,iky,ikz)*conjg(ttwv(ikx,iky,ikz)) )
                 spz(j,14) = spz(j,14) + vzx/wk**2 + vzy/wk**2 + vzz/wk**2
                 spz(j,15) = spz(j,15) + vt*aj/bj
              endif
           enddo
        enddo
     enddo
  endif

  if (potkinwv_flag == 1) then
     deallocate(gewv)
     call mpi_reduce(spz,spztot,nspz15,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,istatus)
  else
     call mpi_reduce(spz,spztot,nspz13,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,istatus)
  endif
  call mpi_reduce(  n,  ntot,  nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,istatus)
     
  if (mype.eq.0) then
     do j=j0,kts-1 
        if (ntot(j).ne.0) then
           if (potkinwv_flag == 1) then
              write(iu,5011) float(j),(spztot(j,i),i=1,15),ntot(j)
           else
              write(iu,5000) float(j),(spztot(j,i),i=1,13),ntot(j)
           endif
        endif
     enddo
     write(iu,*) '           '
     !    write(iu,*) '           '
    call flush(iu)
 endif
 deallocate(spz,spztot)
 
 return
5000 format(1X,F4.0,4X,13(E13.6,1x),6X,I12)
5011 format(1X,F4.0,4X,15(E13.6,1x),6X,I12)
end subroutine spec


subroutine transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k,nuk,nvk,nwk,ispec,iu)

! Calculates transfer spectra.
! ispec = 1: k
! ispec = 2: kh
! ispec = 3: kz

  implicit none 
  include 'mpif.h'

  integer, intent(in) :: ispec,iu
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,nzx,nzy,nzz,ntt,geok,gw1k,gw2k
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: ngeok,ngw1k,ngw2k,nuk,nvk,nwk

  integer :: ikx,iky,ikz,ikza,j,j0,nspz,nspz6
  integer, dimension(0:kts) :: n,ntot 
  real :: kx,ky,kz,k,k2,kh2,wkh,kh2n,kzn
  real :: vzx,vzy,vzz,vtt
  real, dimension(0:kts,6)  :: spz,spztot
  complex :: u,v,w,c1,c2,c3

  nspz  = kts+1
  nspz6 = (kts+1)*6

  if (ispec.eq.1) then
    j0 = 1
  elseif (ispec.eq.2) then
    j0 = 0
  elseif (ispec.eq.3) then
    j0 = 0
  else
    print*,"ispec error"
    stop
  endif
  
  call wtoab(nzx,nzy,nzz,ntt,ngeok,ngw1k,ngw2k,nuk,nvk,nwk)
  call velo(nzx,nzy,nzz,nuk,nvk,nwk)

  N   = 0
  spz = 0.

  do ikz = 1,iktzp
    ikza = mype*iktzp+ikz
    kz = kza(ikza)
    kzn = kz*L3/twopi
    do iky = 1,ikty
      ky  = kya(iky)
      do ikx = 1,iktx
        kx   = kxa(ikx)
        kh2  = kx*kx + ky*ky
        kh2n = kh2 * (L1/twopi)**2
        kh2  = max(1.e-15,kh2)
        wkh  = sqrt(kh2)
        k2   = kx*kx + ky*ky + kz*kz
        k    = sqrt(k2)
        k2   = max(1.e-15,k2)
           
        if (ispec.eq.1) then
          j = int(k*L1/twopi+0.5)
        elseif (ispec.eq.2) then
          j = int(wkh*L1/twopi+0.5)
        elseif (ispec.eq.3) then
          j = int(abs(kz)*L3/twopi+0.5)
        endif

        if (L(ikx,iky,ikz).eq.1)  then
          if (j.lt.j0 .or. j.gt.kts) then
            print*,'transf: screw-up.   k= ',j,kx,ky,kz,L(ikx,iky,ikz)
          endif

          c1 = ky*zz(ikx,iky,ikz) - kz*zy(ikx,iky,ikz)
          c2 = kz*zx(ikx,iky,ikz) - kx*zz(ikx,iky,ikz)
          c3 = kx*zy(ikx,iky,ikz) - ky*zx(ikx,iky,ikz)
          u = zi * c1 / k2
          v = zi * c2 / k2
          w = zi * c3 / k2

          ! k \in R_k
          if (kh2n.gt.1.e-10 .and. abs(kzn).gt.1.e-10) then
            spz(j,1) = spz(j,1) + real(geok(ikx,iky,ikz)*conjg(ngeok(ikx,iky,ikz)))/kh2
            spz(j,2) = spz(j,2) + real(gw1k(ikx,iky,ikz)*conjg(ngw1k(ikx,iky,ikz)))/kh2
            spz(j,2) = spz(j,2) + real(gw2k(ikx,iky,ikz)*conjg(ngw2k(ikx,iky,ikz)))/kh2
          endif

          ! Special cases: i) k_h=0, ii) k_z=0.
          vzx = real(u              *conjg( nuk(ikx,iky,ikz)))
          vzy = real(v              *conjg( nvk(ikx,iky,ikz)))
          vzz = real(w              *conjg( nwk(ikx,iky,ikz)))
          vtt = real(tt(ikx,iky,ikz)*conjg(ntt(ikx,iky,ikz)))

          !  k \in V_k
          if (kh2n.lt.1.e-10.and.abs(kzn).gt.1.e-10) then
            spz(j,1) = spz(j,1) + vtt
            spz(j,2) = spz(j,2) + vzx + vzy            
          endif

          !  k \in B_k
          if (kh2n.gt.1.e-10.and.abs(kzn).lt.1.e-10) then
            spz(j,1) = spz(j,1) + vzx + vzy 
            spz(j,2) = spz(j,2) + vzz + vtt*aj/bj
          endif

          !  Now, compute KE and PE transfer
          spz(j,4) = spz(j,4) + vzx + vzy + vzz
          spz(j,5) = spz(j,5) + vtt*aj/bj

          n(j)   = n(j) + 2
        endif
      enddo
    enddo
  enddo

  spz(:,3) = spz(:,1) + spz(:,2)
  spz(:,6) = spz(:,4) + spz(:,5)

  call mpi_reduce(spz,spztot,nspz6,MPI_REAL,   MPI_SUM,0,MPI_COMM_WORLD,istatus)
  call mpi_reduce(  n,  ntot, nspz,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,istatus)

  if (mype.eq.0) then  
    do j = j0,kts-1
      if (ntot(j).ne.0) then
        write(iu,5000) float(j),(spztot(j,iky),iky=1,6),ntot(j)
      endif
    enddo
    write(iu,*) '     '
    ! write(iu,*) '     '
  endif

  return
5000 format(1X,F4.0,4X,6(E11.4,1x),10X,I6)
end subroutine transf

subroutine do_diagnostics(zx,zy,zz,tt,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,rhtt,uk,vk,wk,nt)
  !! Do all the diagnostics here. Customise based on what you need for analysis
  implicit none
  include 'mpif.h'
  integer, intent(in)  :: nt
  complex, intent(in), dimension(iktx,ikty,iktzp) :: zx,zy,zz,tt,nzxk,nzyk,nzzk,nttk,geok,gw1k,gw2k
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: rhzx,rhzy,rhzz,rhtt,uk,vk,wk
  
  !! Note in the following subroutines (spec,trans,out_eng) zx,zy,zz,tt,nzx,nzy,nzz,ntt,geok,gw1k,gw2k
  !! remain unchanged but 
  call out_eng(zx,zy,zz,tt,uk,vk,wk,geok,gw1k,gw2k,rhzx,rhzy,rhzz,rhtt,nt)
  if (keepSPCH==1) call spec(zx,zy,zz,tt,uk,vk,wk,geok,gw1k,gw2k,rhzx,rhzy,rhzz,rhtt,2,iuSPCH)
  if (keepSPCZ==1) call spec(zx,zy,zz,tt,uk,vk,wk,geok,gw1k,gw2k,rhzx,rhzy,rhzz,rhtt,3,iuSPCZ)
  if (keepSPC==1 ) call spec(zx,zy,zz,tt,uk,vk,wk,geok,gw1k,gw2k,rhzx,rhzy,rhzz,rhtt,1,iuSPC )
  if (keepTRNH==1) call transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,2,iuTRNH)
  if (keepTRNZ==1) call transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,3,iuTRNZ)
  if (keepTRN ==1) call transf(zx,zy,zz,tt,geok,gw1k,gw2k,nzxk,nzyk,nzzk,nttk,rhzx,rhzy,rhzz,uk,vk,wk,1,iuTRN )
  
  return
end subroutine do_diagnostics

 
end module diagnostics

   
