!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module dumps subspaces of real (physcial) space in the file 'realspace.ncf'             !!
!! At the moment is limited to 2D slices in physical(real) but can be extended to 3D subsets    !!
!! It can be customised to store one of the fields below                                        !!
!! vorticity component, temperature, velocity component                                         !!
!! other fields can also be stored                                                              !! 
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module realspacedumps
  use param
  use param_fftw
  use netcdf
  use velvorproj
  use nm_decomp
  implicit none
  
  !! %% PARAMETERs %%
  integer, parameter :: nHor = FLOOR(real(nHend-nHstart+1.0,8)/nskipH)
  integer, parameter :: nVer = FLOOR(real(nVend-nVstart+1.0,8)/nskipV)

  ! the location of slices
  integer, parameter :: ihor = n3/2, iver = n1d/2
  
  ! Here we determine which fields are going to be stored
  ! Vorticity:
  ! horizontal slices
  integer, parameter :: kpZZHtot=1, kpZZHbrt=0
  ! vertical slices
  integer, parameter :: kpZZVtot=1, kpZZVbrt=0
  ! Velocities
  ! horizontal slices
  integer, parameter :: kpUUHtot=1, kpVVHtot=1, kpUUHbrc=0, kpVVHbrc=0
  ! vertical slices
  integer, parameter :: kpUUVtot=1, kpVVVtot=1, kpUUVbrc=0, kpVVVbrc=0
  ! Quadratic Nonlinear terms to derive pressure
  ! horizontal slices
  integer, parameter :: kpQDHtot=0, kpQDHbrc=0
  ! vertical slices
  integer, parameter :: kpQDVtot=0, kpQDVbrc=0  
  
  !! %% VARIABLES %%
  ! NetCDF IDs
  integer, save  :: idslices,idxrs,idyrs,idzrs,idtmrs,tmrid,xxid,yyid,zzid
  integer, save  :: zzhtotid,zzvtotid,zzhbrtid,zzvbrtid
  integer, save  :: uuhtotid,uuvtotid,vvhtotid,vvvtotid
  integer, save  :: uuhbrcid,uuvbrcid,vvhbrcid,vvvbrcid
  integer, save  :: qdhtotid,qdvtotid,qdhbrcid,qdvbrcid
  integer, save  :: iRScount   ! counter when real space fields are dumped
   

  !! Make internal variables and functions private
  PRIVATE        :: kpZZHtot,kpZZHbrt,kpZZVtot,kpZZVbrt
  PRIVATE        :: kpUUHtot,kpVVHtot,kpUUHbrc,kpVVHbrc
  PRIVATE        :: kpUUVtot,kpVVVtot,kpUUVbrc,kpVVVbrc
  PRIVATE        :: kpQDHtot,kpQDHbrc,kpQDVtot,kpQDVbrc
  PRIVATE        :: check_rs
  
CONTAINS 

subroutine check_rs(sstatus)
  implicit none
  include 'mpif.h'  
  integer, intent (in) :: sstatus
  if(sstatus /= nf90_noerr) then
     print *, 'ERROR in realspacedumps.F90'
     print *, trim(nf90_strerror(sstatus))
     call MPI_Abort(MPI_COMM_WORLD)
  end if
end subroutine check_rs
  
subroutine prep_realslice()
  ! creates a netcdf file for slices of fields in real (physical) space.
  implicit none
  include 'mpif.h'
!!$  real, dimension(nHor):: xr,yr
!!$  real, dimension(nVer):: zr
  real    :: xr,yr,zr
  integer :: ix,iy,iz,ix1,iy1,iz1
  
  iRScount = 0
  
  if ((iver.gt.n2).or.(ihor.gt.n3)) then
     print*, '----------- the slice index out of bound -----------------'
  endif
   if (mod(n2dp,nskipH).ne.0)  then
     print*, '----------- cannot handle this skipping: n2dp cannot be divided by nskipH -----------------'
  endif

  istatus =  nf90_create("realspace.ncf",ior(NF90_NETCDF4,NF90_MPIIO), &
       idslices,comm= MPI_COMM_WORLD,info = MPI_INFO_NULL)
  if (istatus.ne.0) print*,'Yo! error in creating realspace.ncf! Darn!'
  
  ! Define the dimensions of space and time
  call check_rs( nf90_def_dim(idslices,"xxrs",nHor,idxrs))
  call check_rs( nf90_def_dim(idslices,"yyrs",nHor,idyrs))
  call check_rs( nf90_def_dim(idslices,"zzrs",nVer,idzrs))
  call check_rs( nf90_def_dim(idslices,"timers",nrsout+1,idtmrs)) ! nsp2dout+ >>1<< for dumping IC
  ! Define the time variable and spatial coordinates
  call check_rs( nf90_def_var(idslices,"vartime",NF90_FLOAT,(/ idtmrs /),tmrid))
  call check_rs( nf90_def_var(idslices,"varxx",NF90_FLOAT,(/ idxrs /),xxid))
  call check_rs( nf90_def_var(idslices,"varyy",NF90_FLOAT,(/ idyrs /),yyid))
  call check_rs( nf90_def_var(idslices,"varzz",NF90_FLOAT,(/ idzrs /),zzid))

  ! Define the variables for the slices that are supposed to be stored
  ! i.e. zx, zy, zz and tt

  if (kpZZHtot.eq.1) then
     istatus=nf90_def_var(idslices,"ZZHtot",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),zzhtotid)
  endif
  if (kpZZVtot.eq.1) then
     istatus=nf90_def_var(idslices,"ZZVtot",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),zzvtotid)
  endif
  if (kpZZHbrt.eq.1) then
     istatus=nf90_def_var(idslices,"ZZHbrt",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),zzhbrtid)
  endif
  if (kpZZVbrt.eq.1) then
     istatus=nf90_def_var(idslices,"ZZVbrt",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),zzvbrtid)
  endif
 
  if (kpUUHtot.eq.1) then
     istatus=nf90_def_var(idslices,"UUHtot",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),uuhtotid)
  endif
  if (kpUUVtot.eq.1) then
     istatus=nf90_def_var(idslices,"UUVtot",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),uuvtotid)
  endif
  if (kpVVHtot.eq.1) then
     istatus=nf90_def_var(idslices,"VVHtot",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),vvhtotid)
  endif
  if (kpVVVtot.eq.1) then
     istatus=nf90_def_var(idslices,"VVVtot",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),vvvtotid)
  endif
  if (kpUUHbrc.eq.1) then
     istatus=nf90_def_var(idslices,"UUHbrc",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),uuhbrcid)
  endif
  if (kpUUVbrc.eq.1) then
     istatus=nf90_def_var(idslices,"UUVbrc",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),uuvbrcid)
  endif
  if (kpVVHbrc.eq.1) then
     istatus=nf90_def_var(idslices,"VVHbrc",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),vvhbrcid)
  endif
  if (kpVVVbrc.eq.1) then
     istatus=nf90_def_var(idslices,"VVVbrc",NF90_FLOAT,(/idyrs,idzrs,idtmrs/),vvvbrcid)
  endif

  if (kpQDHtot.eq.1) then
     istatus=nf90_def_var(idslices,"QDHtot",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),qdhtotid)
  endif
  if (kpQDVtot.eq.1) then
     istatus=nf90_def_var(idslices,"QDVtot",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),qdvtotid)
  endif
  if (kpQDHbrc.eq.1) then
     istatus=nf90_def_var(idslices,"QDHbrc",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),qdhbrcid)
  endif
  if (kpQDVbrc.eq.1) then
     istatus=nf90_def_var(idslices,"QDVbrc",NF90_FLOAT,(/idxrs,idyrs,idtmrs/),qdvbrcid)
  endif
  
  
  call check_rs( nf90_enddef(idslices) )

  ! put values in spatial coordinate variables i.e. varxx and varyy vectors
  if (mype.eq.0) then
     ix1 = 0
     iy1 = 0
     iz1 = 0
     do ix = nHstart, nHend, nskipH
        ix1 = ix1 + 1
        xr = (real(ix)-1.0)/(real(n1)-1.0)*L1-L1/2.0
        istatus =  nf90_put_var(idslices,xxid,xr,start = (/ ix1 /))
        if (istatus.ne.0) print*,'Yo! we fucked up xr! Darn!'
     enddo
     do iy = nHstart, nHend, nskipH
        iy1 = iy1 + 1
        yr = (real(iy)-1.0)/(real(n2d)-1.0)*L2-L2/2.0
        call check_rs( nf90_put_var(idslices,yyid,yr,start = (/ iy1 /)) )
     enddo
     do iz = nVstart, nVend, nskipV
        iz1 = iz1 + 1
        zr = (real(iz)-1.0)/(real(n3d)-1.0)*L3-L3/2.0
        call check_rs( nf90_put_var(idslices,zzid,zr,start = (/ iz1 /)) )
     enddo
  endif

  return
end subroutine prep_realslice

subroutine dump_realspace(zxk,zyk,zzk,ttk,uk,vk,wk,zxr,zyr,zzr,ttr,ur,vr,wr)
  
  implicit none 
  include 'mpif.h'
  real,    intent(inout), dimension(n1d,n3d,n2dp)    :: zxr,zyr,zzr,ttr,ur,vr,wr
  complex, intent(inout), dimension(iktx,ikty,iktzp) :: zxk,zyk,zzk,ttk,uk,vk,wk
  complex, dimension(iktx,ikty,iktzp)   :: zxaux,zyaux,zzaux,ttaux,qdk
  integer :: ikz,ikza,ikx,iky
  real    :: kz,kx,ky

  iRScount = iRScount + 1
  if (mype.eq.0) istatus = nf90_put_var(idslices,tmrid, time, start = (/ iRScount /))

  zzaux = zzk
  zxaux = zxk
  zyaux = zyk
  ttaux = ttk
  
  call fftwkr(plan3_zznk_zznr,zzk,zzr)
  if (kpZZHtot.eq.1) call dump_horslice(zzr,zzhtotid)
  if (kpZZVtot.eq.1) call dump_verslice(zzr,zzvtotid)
  call fftwrk(plan3_zznr_zznk,zzr,zzk) 

  ! First store the total velocities
  call velo(zxk,zyk,zzk,uk,vk,wk)
  call fftwkr(plan3_uk_ur,uk,ur)
  if (kpUUHtot.eq.1) call dump_horslice(ur,uuhtotid)
  if (kpUUVtot.eq.1) call dump_verslice(ur,uuvtotid)
  zxr = ur
  call fftwrk(plan3_ur_uk,ur,uk)
  
  call fftwkr(plan3_vk_vr,vk,vr)
  if (kpVVHtot.eq.1) call dump_horslice(vr,vvhtotid)
  if (kpVVVtot.eq.1) call dump_verslice(vr,vvvtotid)
  zyr = vr
  call fftwrk(plan3_vr_vk,vr,vk)

  ! Derive the quadratic nonlinear terms for calculating pressure (using total vel.)
  qdk = cmplx(0.,0.)
  ttr = zxr*zxr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-kx*kx*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  qdk = qdk + ttk

  ttr = zxr*zyr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-2*kx*ky*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo
  qdk = qdk + ttk

  ttr = zyr*zyr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx  
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-ky*ky*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  qdk = qdk + ttk
  ttk = qdk
  call fftwkr(plan3_ttnk_ttnr,ttk,ttr)
  if (kpQDHtot.eq.1) call dump_horslice(ttr,qdhtotid)
  if (kpQDVtot.eq.1) call dump_verslice(ttr,qdvtotid)
  
  zxk = cmplx(0.,0.)
  zyk = cmplx(0.,0.)
  

  ! Derive the baroclinic velocities and store them
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              uk(ikx,iky,ikz)=cmplx(0.,0.)
              vk(ikx,iky,ikz)=cmplx(0.,0.)
           else
              zzk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  call fftwkr(plan3_zznk_zznr,zzk,zzr)
  if (kpZZHbrt.eq.1) call dump_horslice(zzr,zzhbrtid)
  if (kpZZVbrt.eq.1) call dump_verslice(zzr,zzvbrtid)
  call fftwrk(plan3_zznr_zznk,zzr,zzk)
  
  call fftwkr(plan3_uk_ur,uk,ur)
  if (kpUUHbrc.eq.1) call dump_horslice(ur,uuhbrcid)
  if (kpUUVbrc.eq.1) call dump_verslice(ur,uuvbrcid)
  call fftwrk(plan3_ur_uk,ur,uk) 

  call fftwkr(plan3_vk_vr,vk,vr)
  if (kpVVHbrc.eq.1) call dump_horslice(vr,vvhbrcid)
  if (kpVVVbrc.eq.1) call dump_verslice(vr,vvvbrcid)
  call fftwrk(plan3_vr_vk,vr,vk)

  ! Derive the quadratic nonlinear terms for calculating pressure (using barotropic vel.)
  
  zxk = uk
  zyk = vk
  call fftwkr(plan3_zxnk_zxnr,zxk,zxr)
  call fftwkr(plan3_zynk_zynr,zyk,zyr)
  
  qdk = cmplx(0.,0.)
  ttr = zxr*zxr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-kx*kx*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  qdk = qdk + ttk

  ttr = zxr*zyr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-2*kx*ky*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  qdk = qdk + ttk

  ttr = zyr*zyr
  call fftwrk(plan3_ttnr_ttnk,ttr,ttk)
  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           if (abs(kz).lt.1e-4) then
              ttk(ikx,iky,ikz)=-ky*ky*ttk(ikx,iky,ikz)
           else
              ttk(ikx,iky,ikz)=cmplx(0.,0.)
           endif
        enddo
     enddo
  enddo

  qdk = qdk + ttk
  ttk = qdk
  call fftwkr(plan3_ttnk_ttnr,ttk,ttr)
  if (kpQDHbrc.eq.1) call dump_horslice(ttr,qdhbrcid)
  if (kpQDVbrc.eq.1) call dump_verslice(ttr,qdvbrcid)
  

  zzk = zzaux
  zxk = zxaux  
  zyk = zyaux  
  ttk = ttaux
  
  return
end subroutine dump_realspace

subroutine dump_horslice(fr,idhvar)
  implicit none
  include 'mpif.h'
  
  real, intent(in), dimension(n1d,n3d,n2dp)  :: fr
  real,             dimension(n1d,n2dp)      :: frslice,junk
  real,             dimension(n1d,n2d)       :: horvars
  integer, intent(in)   :: idhvar
  integer, dimension(3) :: nccount,ncstart
  integer               :: k,nsends,nbuf,iproc,istart
  integer               :: status(MPI_STATUS_SIZE)

  nsends = npe - 1
  do k=1,n2dp
     frslice(:,k) = fr(:,ihor,k)
  enddo
  nbuf = n1d*n2dp
  if (mype.gt.0) then
     call mpi_send(frslice,nbuf,MPI_REAL,0,137,MPI_COMM_WORLD,istatus)
  endif

  if (mype.eq.0) then
     horvars(1:n1d,1:n2dp) = frslice   
     do iproc=1,nsends
        call mpi_recv(junk,nbuf,MPI_REAL,MPI_ANY_SOURCE,137,MPI_COMM_WORLD,status,istatus)
        istart=status(MPI_SOURCE)*n2dp
        do k=1,n2dp
           horvars(1:n1d,istart+k)=junk(1:n1d,k)
        enddo
     enddo
     ! print*, 'horvars(384,384) = ', horvars(384,384)
     call check_rs( nf90_put_var(idslices,idhvar, &
          horvars(nHstart+1:nHend+1:nskipH,nHstart+1:nHend+1:nskipH), &
          start= (/ 1, 1, iRScount/), count= (/ nHor, nHor, 1/)))
  endif
  
end subroutine dump_horslice

subroutine dump_verslice(fr,idvvar)
  implicit none
  include 'mpif.h'
  
  real, intent(in), dimension(n1d,n3d,n2dp) :: fr
  integer, intent(in)   :: idvvar
  real,             dimension(nHor,nVer)    :: vervars
  integer   :: wantedmype,ix,iz,ix1,iz1

  wantedmype = floor((iver-1.0)/n2dp)
  
  if (wantedmype.eq.mype) then
     iz1=0
     do iz = nVstart, nVend, nskipV
        iz1=iz1+1
        ix1=0
        do ix = nHstart+1, nHend+1, nskipH
           ix1=ix1+1
           vervars(ix1,iz1) = fr(ix,iz,iver-wantedmype*n2dp)
        end do
     end do
     call check_rs(nf90_put_var(idslices, idvvar, vervars, &
          start= (/ 1, 1, iRScount/),count=(/ nHor, nVer, 1/)))
  endif
  
end subroutine dump_verslice

subroutine close_ncf_realspace()
  ! close creates a netcdf file for real space dumping
  implicit none
  include 'mpif.h'
  call check_rs( nf90_close(idslices) )
  return
end subroutine close_ncf_realspace

end module realspacedumps
 
