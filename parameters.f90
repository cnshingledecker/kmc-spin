MODULE parameters
  USE typedefs
  SAVE
  !******************************************************************************
  ! Main Parameters
  !******************************************************************************
  INTEGER         , PARAMETER :: L        = 100     ! Lattice size
  INTEGER         , PARAMETER :: EMAX     = 1000     ! Max number of electrons allowed in model at any time
  DOUBLE PRECISION, PARAMETER :: MAXTIME  = 10.0d0 ! The maximum model runtime
  DOUBLE PRECISION, PARAMETER :: T_KELVIN = 300.0  ! The temperature of the lattice (K)
  DOUBLE PRECISION, PARAMETER :: VOLTAGE  = 0.0    ! Applied Voltage
  DOUBLE PRECISION, PARAMETER :: CHFI     = 1.0    ! HFI coupling constant
  DOUBLE PRECISION, PARAMETER :: CJ       = 0.0    ! Exchange coupling constant
  DOUBLE PRECISION, PARAMETER :: CSOC     = 0.0    ! SOC coupling constant
  DOUBLE PRECISION, PARAMETER :: KHOP     = 30     ! 1.0d0/(DBLE(NSTEPS)*DT)
  DOUBLE PRECISION, PARAMETER :: FIXEDTAU = (1.0d0/KHOP) ! (MAXTIME/REAL(NSTEPS))
  DOUBLE PRECISION, PARAMETER :: DT       = FIXEDTAU/100.0 ! Time step (Delta t)
  DOUBLE PRECISION            :: NEXTHOP  = 0.0d0  ! The time of the next hop
  DOUBLE PRECISION            :: TAUFAC   = 1 
  INTEGER         , PARAMETER :: NSTEPS   = MAXTIME/FIXEDTAU !INT(1.0d0/(KHOP*DT))/6 ! Number of time steps between hops

  !******************************************************************************
  ! Physical Conditions
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: SIGMA_DOS  = 1.0      ! Width of DOS
  DOUBLE PRECISION, PARAMETER :: EPS_R      = 3.0      ! \epsilon_r
  DOUBLE PRECISION, PARAMETER :: A_NN       = 1.6E-9   ! Distance to nearest neighbor
  DOUBLE PRECISION, PARAMETER :: ALPHA      = 10./A_NN ! See Cottaar thesis
  INTEGER         , PARAMETER :: N_IMAGES   = 100      ! See Cottaar thesis
  INTEGER         , PARAMETER :: N_DISKS    = 100000   ! See Cottaar thesis
  DOUBLE PRECISION, PARAMETER :: R_c        = 0.0d0 !10.*A_NN ! Critical radius for explicit Coulomb calculations

  !******************************************************************************
  ! Equivalent of Common Block: TIME
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: DT2  = DT/2.d0
  DOUBLE PRECISION, PARAMETER :: DT3  = DT/3.d0
  DOUBLE PRECISION, PARAMETER :: DT6  = DT/6.d0
  DOUBLE PRECISION, PARAMETER :: DT24 = DT/24.d0

  !******************************************************************************
  ! Spin/MD Parameters (from md.inp and test.inp)
  !******************************************************************************
  INTEGER         , PARAMETER :: NNEIB     = 6 ! Number of nearest neighbors
  DOUBLE PRECISION, PARAMETER :: DISTNN    = 1.0 !Nearest neighbor distance
  DOUBLE PRECISION, PARAMETER :: MASS      = 1.0 !Mass
  DOUBLE PRECISION, PARAMETER :: MBSCALE   = DSQRT(T_KELVIN/MASS)

  !******************************************************************************
  ! Spin/MD Parameters (from code)
  !******************************************************************************
  INTEGER                     :: NSPIN     = 0 ! Number of electrons in model at any point
  DOUBLE PRECISION, PARAMETER :: EPS       = 1.0d0
  DOUBLE PRECISION, PARAMETER :: DLT       = 1.0d0
  DOUBLE PRECISION, PARAMETER :: SIG       = 1.0d0
  DOUBLE PRECISION, PARAMETER :: CELL      = DBLE(L)*DISTNN
  DOUBLE PRECISION, PARAMETER :: RCUT      = 0.5*CELL
  DOUBLE PRECISION, PARAMETER :: RLIST     = 1.05*RCUT
  DOUBLE PRECISION, PARAMETER :: RCUT2     = RCUT*RCUT
  DOUBLE PRECISION, PARAMETER :: RLIST2    = RLIST*RLIST
  DOUBLE PRECISION, PARAMETER :: DLTR      = RCUT/250.d0

  !******************************************************************************
  ! Matrix/Crystal Structure Parameters
  !******************************************************************************
  TYPE (node), ALLOCATABLE, TARGET  :: MATRIX(:,:,:)  ! Ice-mantle matrix
  INTEGER, PARAMETER                :: NEDGE     =  L !Edge elements
  INTEGER, PARAMETER                :: NTHICK    =  L !Thickness elements
  INTEGER                           :: NELEMS    =  0 ! Number of elements in lattice
  INTEGER, PARAMETER                :: DIMENS(3) = (/ NTHICK, NEDGE, NEDGE /)

  !******************************************************************************
  ! Kinetic Parameters
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: NU_0       = 2.6E11   ! Trial frequency, for the rates, in 1/s

  !******************************************************************************
  ! Constants
  !******************************************************************************
  DOUBLE PRECISION, PARAMETER :: PI        = 4.D0*DATAN(1.D0) ! Pi
  DOUBLE PRECISION, PARAMETER :: EBASE     = EXP(1.D0)        ! Natural base
  DOUBLE PRECISION, PARAMETER :: EPS_0     = 8.854E-12        ! C/(V*m)
  DOUBLE PRECISION, PARAMETER :: K_B       = 8.617332E-5      ! eV/K
  DOUBLE PRECISION, PARAMETER :: ECHARGE   = 1.602176E-19     ! Coulombs
  DOUBLE PRECISION, PARAMETER :: VC2EV     = 6.2415091E18     ! Volt Coulombs to eV

  !******************************************************************************
  ! Counters
  !******************************************************************************
  DOUBLE PRECISION            :: TIME            ! Global time counter
  INTEGER         , PARAMETER :: TIME_FREQ  = 1  !1000000
  INTEGER         , PARAMETER :: E_ADD      = 5
  INTEGER                     :: ECOUNT     = 0

  !******************************************************************************
  ! I/O Unit Numbers
  !******************************************************************************
  INTEGER         , PARAMETER :: LAYERUNIT  = 1001
  INTEGER         , PARAMETER :: MATRIXUNIT = 1002

  !******************************************************************************
  ! Lattice locations for array elements
  !******************************************************************************
  INTEGER         , ALLOCATABLE :: SPINLOCS(:,:)

  !******************************************************************************
  ! Equivalent of Common Block: Rxyz
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: X(:)
  DOUBLE PRECISION, ALLOCATABLE :: Y(:)
  DOUBLE PRECISION, ALLOCATABLE :: Z(:)

  !******************************************************************************
  ! Equivalent of Common Block: Vxyz
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: VX(:)
  DOUBLE PRECISION, ALLOCATABLE :: VY(:)
  DOUBLE PRECISION, ALLOCATABLE :: VZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: Fxyz
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: FX(:)
  DOUBLE PRECISION, ALLOCATABLE :: FY(:)
  DOUBLE PRECISION, ALLOCATABLE :: FZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: S
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SX(:)
  DOUBLE PRECISION, ALLOCATABLE :: SY(:)
  DOUBLE PRECISION, ALLOCATABLE :: SZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: S0
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SX0(:)
  DOUBLE PRECISION, ALLOCATABLE :: SY0(:)
  DOUBLE PRECISION, ALLOCATABLE :: SZ0(:)

  !******************************************************************************
  ! Equivalent of Common Block: S1
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SX1(:)
  DOUBLE PRECISION, ALLOCATABLE :: SY1(:)
  DOUBLE PRECISION, ALLOCATABLE :: SZ1(:)

  !******************************************************************************
  ! Equivalent of Common Block: S2
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SX2(:)
  DOUBLE PRECISION, ALLOCATABLE :: SY2(:)
  DOUBLE PRECISION, ALLOCATABLE :: SZ2(:)

  !******************************************************************************
  ! Equivalent of Common Block: S3
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SX3(:)
  DOUBLE PRECISION, ALLOCATABLE :: SY3(:)
  DOUBLE PRECISION, ALLOCATABLE :: SZ3(:)

  !******************************************************************************
  ! Equivalent of Common Block: WS
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: WSX(:)
  DOUBLE PRECISION, ALLOCATABLE :: WSY(:)
  DOUBLE PRECISION, ALLOCATABLE :: WSZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: WK
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: WKX(:)
  DOUBLE PRECISION, ALLOCATABLE :: WKY(:)
  DOUBLE PRECISION, ALLOCATABLE :: WKZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: B
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: BX(:)
  DOUBLE PRECISION, ALLOCATABLE :: BY(:)
  DOUBLE PRECISION, ALLOCATABLE :: BZ(:)

  !******************************************************************************
  ! Equivalent of Common Block: SOCXYZ
  !******************************************************************************
  DOUBLE PRECISION, ALLOCATABLE :: SOCX(:)
  DOUBLE PRECISION, ALLOCATABLE :: SOCY(:)
  DOUBLE PRECISION, ALLOCATABLE :: SOCZ(:)

  !******************************************************************************
  ! Switches
  !******************************************************************************
  !******************************************************************************
  ! Select NO_HOP to populate every lattice site with an electron
  !******************************************************************************
  LOGICAL         , PARAMETER :: NO_HOP        = .FALSE.

  !******************************************************************************
  ! Select DEBUG to enable extra print statements for correctness checking
  !******************************************************************************
  LOGICAL         , PARAMETER :: DEBUG         = .FALSE.

  !******************************************************************************
  ! Select VERBOSE for extra print statements (fewer print statements than debug)
  !******************************************************************************
  LOGICAL         , PARAMETER :: VERBOSE       = .FALSE.

  !******************************************************************************
  ! Select SOC_ON to turn spin-orbit-coupling on
  !******************************************************************************
  LOGICAL         , PARAMETER :: SOC_ON        = .FALSE.

  !******************************************************************************
  ! Select MD_ON to enable extra MD data structures, e.g.  vi, fi, etc.
  !******************************************************************************
  LOGICAL         , PARAMETER :: MD_ON         = .FALSE.

  !******************************************************************************
  ! Select INFINITE_HOPS to make x-coord periodic,
  ! i.e. spin information isn't lost to the electrodes
  !******************************************************************************
  LOGICAL         , PARAMETER :: INFINITE_HOPS = .TRUE.

  !******************************************************************************
  ! Select INITIAL_ELECS to populate the lattice with EMAX (see aboce) electrons
  ! at random initial positions in the latttice before the start of the main loop
  !******************************************************************************
  LOGICAL         , PARAMETER :: INITIAL_ELECS = .FALSE.

  !******************************************************************************
  ! Select FIXEDSTEP to enable constant time between electron hops. This allows for
  ! better comparison with analytic results
  !******************************************************************************
  LOGICAL         , PARAMETER :: FIXEDSTEP     = .TRUE.
   
END MODULE parameters
