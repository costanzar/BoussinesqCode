module param
  use, intrinsic :: iso_c_binding
  implicit none

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@@@@@@@@@@@@@@ PARAMETERS AND GLOBAL VARIABLES @@@@@@@@@@@@@@@@@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


! =================================================================
! ================ SIMULATION RELATED PARAMETERS ==================
! =================================================================

! -----------------------------------------------------------------
! Set Model Resolution
  integer(C_INTPTR_T), parameter :: n1=2048, n2=2048, n3=256
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Number of Processors (must divide n2 and iktz)
  integer, save         :: npe = 256
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Math constants
  real, parameter            :: twopi = 4.*asin(1.)
  real, parameter            :: sqrt2 = sqrt(2.)
  complex, parameter :: zi    = cmplx(0.,1.)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Time
  real, save         :: time                                       ! global variable showing time
  real,    parameter :: delt  =  1.7                              ! timestep
  real,    parameter :: tstop =  1*24*3600                         ! length of integration
  integer, parameter :: nstop  = int(tstop/delt)                   ! number of timesteps

! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Stratification and Rotation
  real, parameter :: PrRa =   100                                  ! Prandtle Ratio  N/f
  real, parameter :: aj   =   1e-2                                 ! thermal expansivity, N^2 = aj*bj
  real, parameter :: bj   =   aj                                   ! background theta gradient
  real, parameter :: bf2  =   aj*bj                                ! Brunt-Vaisalla freq squared
  real, parameter :: bf   =   sqrt(bf2)                            ! Brunt-Vaisalla frequency
  real, parameter :: cor  =   1.1e-4                               ! Coriolis parameter
  real, parameter :: cor2 =   cor*cor                              ! Coriolis parameter squared
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Set Domain Size
  real, parameter :: L1=10e6, L2=10e6
  real, parameter :: L3=2e4                                ! makes dx=dz
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Physical and Fourier Array Dimensions
  integer(C_INTPTR_T),parameter :: n1d = n1+2, n2d=n2, n3d=n3      ! size of physical arrays (padded for in-place transforms)
  integer(C_INTPTR_T),save      :: n2dp, n2p                       ! size of local arrays with mpi, n2dp=n2d/npe, n2p=n2/npe
  integer(C_INTPTR_T),parameter :: ktx = n1/2, kty =n2/2, ktz =n3/2! max integer wavenumber
  integer(C_INTPTR_T),parameter :: iktx= ktx+1, ikty=n2, iktz=n3   ! no. of wavenumbers
  integer(C_INTPTR_T),save      :: iktzp                           ! no. of local wavenumbers(mpi), iktzp=iktz/npe
  integer(C_INTPTR_T),parameter :: kts = n1                        ! for spectra; should be max(ktx,ktz)
  real, parameter               :: ktrunc_x=twopi/L1*float(n1)/3.  ! dimensional trunc wavenumber (x) not index
  real, parameter               :: ktrunc_y=twopi/L2*float(n2)/3.  ! dimensional trunc wavenumber (y)
  real, parameter               :: ktrunc_z=twopi/L3*float(n3)/3.  ! dimensional trunc wavenumber (z)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Dissipation
  integer, parameter :: ilap = 4, ilap2 = 2*ilap                   ! hyperviscosity order (ilap=1 regular visc)

! -----------------------------------------------------------------
  real, parameter    :: vischtmp = 5.e-5 * (10./ktrunc_x ) **2.
  real, parameter    :: viscztmp = 5.e-5 * (10./ktrunc_z ) **2.
  real, parameter    :: visch = 6e25  !vischtmp * ktrunc_x **(2-2*ilap) ! viscosity coeff horizontal
  real, parameter    :: viscz = 1e5    !!viscztmp * ktrunc_z **(2-2*ilap) ! viscosity coeff vertical
! -----------------------------------------------------------------
! Linear Damping
  real, parameter    :: ek       =  1./(24.*3600.)*delt            ! kh=0 damping time scale
  real, parameter    :: kdamping =  3.                             ! damping on k <= kdamping
! -----------------------------------------------------------------
! Forcing
  integer, parameter :: forcing_flag = 1                           ! = 1 force (see forcing.F90)
! -----------------------------------------------------------------


! -----------------------------------------------------------------
! Wavenumbers
  integer, dimension(:,:,:), allocatable,save :: L                ! Mask for wavenumber truncation
  real, dimension(:), allocatable, save       :: kxa, kya, kza    ! Wavenumbers
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! MPI Stuff
  integer, save :: mype
  integer(C_INTPTR_T), save :: locz, loczstart                     ! size and position of block in z
  integer(C_INTPTR_T), save :: alloc_local                         ! local data size for malloc
  integer, save :: istatus                                         ! used for mpi initiation
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Initial Condition Parameters
  ! icmode = 0 ---> generate IC only based on a descriptions (see init_condition.F90)
  ! icmode = 1 ---> read IC from a netCDF file
  ! icmode = 2 ---> a superpostion of the netCDF file and describing functions
  integer, parameter :: icmode =   1
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! Misc Stuff
  integer, parameter :: iuRESULT = 11                              ! ID for outputting the results (run.list)
  integer, parameter :: truncmode = 0                              ! 1 spherical truncation; 0 cylindrical
  real,    parameter :: fftnorm = float(n1*n2*n3)                  ! used to normal fft/fftinverse
! -----------------------------------------------------------------



! =================================================================
! =============== DIAGNOSTICS RELATED PARAMETERS ==================
! =================================================================
! >> the most frequently-changed parameters of modules are below
! >> the other parameters should be changed inside each modules

! ±±±±±±±±±±±±±±±±±±±± Input and Output Files ±±±±±±±±±±±±±±±±±±±±±
! -----------------------------------------------------------------
  integer, parameter :: ncwrite = 1                                ! flag for dumping field for other sim's restart
                                                                   ! Note: just the last point in time is dumped

! ±±±±±±±±±±±±±±±±±±±±±±± diagnostics.F90 ±±±±±±±±±±±±±±±±±±±±±±±±±
! >>>  output files: spc?.dat, eng.dat, trn?.dat, run.list and etc
! -----------------------------------------------------------------
! frequency of outputing diagnostics
  integer, parameter :: ndiagout = 48                               ! number of times data dumped
  integer, parameter :: ndiagevery = floor(nstop/float(ndiagout))  ! do diag. every ndiagevery timesteps
! -----------------------------------------------------------------

! -----------------------------------------------------------------

! ±±±±±±±±±±±±±±±± slices of real(physical) space ±±±±±±±±±±±±±±±±±
! >>>  parameters in realspacedumps.F90
! -----------------------------------------------------------------
  integer, parameter :: realspace_flag = 1                         ! = 0 nothing dumped & realspacedumps.F90 not used   
  integer, parameter :: nrsout = 40                                ! number of time slices dumped  
  integer, parameter :: nrsevery=floor(nstop/float(nrsout))        ! dump slices every nrsevery timesteps
  ! starting and end points in slices
  integer, parameter :: nHstart=1, nHend=n1
  integer, parameter :: nVstart=1, nVend=n3
  ! number points to skip
  integer, parameter :: nskipH =  2 ! skip (horizontal) physical points every nskipH (nskip=1 keep all) 
  integer, parameter :: nskipV =  2 ! vertical
! -----------------------------------------------------------------


end module param
