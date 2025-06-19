!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
!! This module contains Normal Mode Decomposition as described in Bartello, 1995, J. Atmos. Sci.  !!
!! The subroutines in this module: wtoab, atowb                                                   !!
!! Modules used : param.F90 and velvorproj.F90                                                    !!
!! Required for: ncf2Dspc.F90 boussinesq.F90 (main file)                                          !!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

module misc
  use param
  use param_fftw
  implicit none

CONTAINS

subroutine printparameters()
  if (mype.eq.0) then
     write(iuRESULT,*) '                '
     write(iuRESULT,*) '                '
     write(iuRESULT,*) 'Parameters: ----------------------------------------------'
     write(iuRESULT,*) '              N1,N2,N3 =  ', n1,n2,n3
     write(iuRESULT,*) '              L1,L2,L3 =  ', L1,L2,L3
     write(iuRESULT,*) '              dx,dy,dz =  ', L1/n1,L2/n2,L3/n3
     write(iuRESULT,*) '              IKTX,Y,Z =  ', iktx,ikty,iktz
     write(iuRESULT,*) '                     KT = ', ktrunc_x,ktrunc_y,ktrunc_z
     write(iuRESULT,*) 'Order of Laplacian diss.=  ', ilap 
     write(iuRESULT,*) '                  VISCH = ', visch
     write(iuRESULT,*) '                  VISCZ = ', viscz
     write(iuRESULT,*) '              tau_VISCH = ', (visch*ktrunc_x**ilap2)**(-1)
     write(iuRESULT,*) '              tau_VISCZ = ', (viscz*ktrunc_z**ilap2)**(-1)
     write(iuRESULT,*) '               Timestep = ', delt
     write(iuRESULT,*) '     Integration length = ', nstop*delt,' = ', nstop,' DT.'
     write(iuRESULT,*) '       Output frequency = ', ndiagevery*delt, ' = ', ndiagevery,' dt.'
     write(iuRESULT,*) '    '
     write(iuRESULT,*) '    Thermal expansivity = ', aj
     write(iuRESULT,*) '    Vertical T gradient = ', bj
     write(iuRESULT,*) '    Brunt-Vaisala freq. = ', sqrt(bf2)
     write(iuRESULT,*) '     Coriolis parameter = ', cor
     write(iuRESULT,*) '                '
     write(iuRESULT,*) '                '
     call flush(iuRESULT)
endif
end subroutine printparameters

!!!!!!! unused routine for now kept here:

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Misc subroutines for fast manifold (move to a module later)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  pvcomp(zxk,zyk,zzk,ttk,uk,vk,wk,zxr,zyr,zzr,ttr,ur,vr,wr)
  ! calculating each term in Pv anamoly in real space and derivin the r.m.s
  ! note here n?k are derivatives of buoyance in x,y and z (not nonlinear terms)
  ! Hence make sure pvcomp is used when changing n?k does not affect the code
  ! like after at the end of all diagnostics such as spec and transf ...
  implicit none 
  
  real,    intent(inout), dimension(n1d,n3d,n2dp)      :: zxr,zyr,zzr,ttr
  real,    intent(inout), dimension(n1d,n3d,n2dp)      :: ur,vr,wr
  complex, intent(inout), dimension(iktx,ikty,iktzp)   :: zxk,zyk,zzk,ttk
  complex, intent(inout), dimension(iktx,ikty,iktzp)   :: uk,vk,wk
  complex,                dimension(iktx,ikty,iktzp)   :: axk,ayk,azk,bxk,byk,bzk

  integer :: ikx,iky,ikz,ikza
  integer :: ii,jj,kk
  real :: kx,ky,kz,term1a,term1,term2,term3,term4,allterms,tmp,maxzz
  complex :: bx,by
  external :: fftwrk,fftwkr
  
  axk = zxk
  ayk = zyk
  azk = zzk
  bxk = uk
  byk = vk
  bzk = wk

  do ikz=1,iktzp
     ikza = mype*iktzp+ikz
     kz = kza(ikza)
     do iky=1,ikty
        ky = kya(iky)
        do ikx=1,iktx
           kx = kxa(ikx)
           uk(ikx,iky,ikz) =  ZI * kx * ttk(ikx,iky,ikz) * bf
           vk(ikx,iky,ikz) =  ZI * ky * ttk(ikx,iky,ikz) * bf
           wk(ikx,iky,ikz) =  ZI * kz * ttk(ikx,iky,ikz) * bf
        enddo
     enddo
  enddo
 
  call fftwkr(plan3_uk_ur,uk,ur) 
  call fftwkr(plan3_vk_vr,vk,vr) 
  call fftwkr(plan3_wk_wr,wk,wr)
  call fftwkr(plan3_zxnk_zxnr,zxk,zxr) 
  call fftwkr(plan3_zynk_zynr,zyk,zyr) 
  call fftwkr(plan3_zznk_zznr,zzk,zzr)


!!$  term1=0
!!$  do ii= 1,n1d
!!$     do jj= 1,n3d
!!$        do kk= 1,n2dp
!!$           term1 = term1 + zzr(ii,jj,kk)*zzr(ii,jj,kk)
!!$        enddo
!!$     enddo
!!$  enddo

  term1 = sum ((zzr/cor)*(zzr/cor))
  term2 = sum ((wr/bf2)*(wr/bf2))
  term3 = sum (((wr/bf2)*(zzr/cor))*((wr/bf2)*(zzr/cor)))
  term4 = sum (((ur/bf2)*(zxr/cor)+(vr/bf2)*(zyr/cor))* &
       ((ur/bf2)*(zxr/cor)+(vr/bf2)*(zyr/cor)))
 
  
  allterms=sum((zzr/cor+wr/bf2+(wr/bf2)*(zzr/cor)+(ur/bf2)*(zxr/cor)+(vr/bf2)*(zyr/cor)) &
       *(zzr/cor+wr/bf2+(wr/bf2)*(zzr/cor)+(ur/bf2)*(zxr/cor)+(vr/bf2)*(zyr/cor)))
  
  call mpi_reduce(term1,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus)
  term1 = sqrt(tmp/float(n1*n2*n3))
  call mpi_reduce(term2,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus)
  term2 = sqrt(tmp/float(n1*n2*n3))
  call mpi_reduce(term3,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus)
  term3 = sqrt(tmp/float(n1*n2*n3))
  call mpi_reduce(term4,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus)
  term4 = sqrt(tmp/float(n1*n2*n3))
  call mpi_reduce(allterms,tmp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,istatus)
  allterms = sqrt(tmp/float(n1*n2*n3))
     
  
  if (mype.eq.0) then ! prep for output
     write(81,5045) time,term1,term2,term3,term4,allterms
     call flush(81)
  
  endif
  
!!$  call fftwrk(plan3_zxnr_zxnk,zxr,zxk) 
!!$  call fftwrk(plan3_zynr_zynk,zyr,zyk) 
!!$  call fftwrk(plan3_zznr_zznk,zzr,zzk)

!!$   if (mype.eq.0) then
!!$     print*,' max of diff = ', maxval(real((zzk-azk)*conjg(zzk-azk)))
!!$   endif
  
  zxk = axk 
  zyk = ayk
  zzk = azk
  uk = bxk
  vk = byk
  wk = bzk

  return

5045  format(1x,f12.2,2x,5(e11.4,1x)) 
end subroutine pvcomp



end module misc
