MODULE subroutines
  USE parameters
  USE bsimple
  USE typedefs
  USE mc_toolbox
  IMPLICIT NONE
CONTAINS
  FUNCTION charge_density ()
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  layer averaged charge density ρ, to be used as a source
    !  term for a calculation of the layer averaged Coulombic
    !  energy, determined using a Poisson equation.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160921      C. Shingledecker     Original code
    !
    !! CHARGE_DENSITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NTHICK) :: charge_density
    INTEGER :: ii, jj, kk

    charge_density = 0
    DO kk = 2, NTHICK-1
       DO jj = 1, NEDGE
          DO ii = 1, NEDGE
             IF ( matrix(ii,jj,kk)%secondary_species .NE. 0 ) THEN
                charge_density(kk-1) = charge_density(kk-1) + 1
             END IF
          END DO
       END DO
       charge_density(kk-1) = charge_density(kk-1)*ECHARGE*(NEDGE**2)
    END DO
  END FUNCTION charge_density

  SUBROUTINE coulomb_sphere (i_x, i_y, i_z, rho_charge)
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate the
    !  the Coulombic energy of a sphere of radius R_c around
    !  some point at i_x, i_y, i_z
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160920      C. Shingledecker     Original code
    !
    !! COULOMB_SPHERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER :: i_x, i_y, i_z ! Coordinates on lattice to calculate energy of
    INTEGER :: j_x, j_y, j_z ! Coordinates on current site
    INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
    INTEGER :: xx, yy, zz
    DOUBLE PRECISION :: R_ij ! Distance between two points on the lattice
    DOUBLE PRECISION :: sr_e, disk_e
    DOUBLE PRECISION, DIMENSION(NTHICK) :: rho_charge
    INTEGER, DIMENSION(6) :: sphere_limits


    ! In the case where the crystal lattice isn't deformed with MD
    sphere_limits = rc_int_dist(i_x, i_y, i_z)
    !    IF ( DEBUG .EQV. .TRUE. ) THEN
    !      PRINT *, "i-coordinates are:",i_x,i_y,i_z
    !      PRINT *, "Sphere limits are:",sphere_limits
    !    END IF

    z_max = sphere_limits(1)
    z_min = sphere_limits(2)
    y_max = sphere_limits(3)
    y_min = sphere_limits(4)
    x_max = sphere_limits(5)
    x_min = sphere_limits(6)

    sr_e = 0.0
    disk_e = 0.0
    DO xx=x_min,x_max
       j_x = i_x + xx
       DO yy=y_min,y_max
          j_y = i_y + yy
          DO zz=z_min,z_max
             j_z = i_z + zz
             ! IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Now in coulomb sphere at:",j_x, j_y, j_z
             IF ( matrix(j_x, j_y, j_z)%secondary_species .EQ. 0 ) THEN
                CONTINUE
             ElSE
                ! Again, in the case where all of the crystal lattice sites are undisturbed by MD
                R_ij = ABS(SQRT((A_NN**2)*(((j_x - i_x)**2) + ((j_y - i_y)**2) + ((j_z - i_z)**2))))
                IF ( R_ij .GT. R_c) THEN
                   CONTINUE
                ELSE
                   ! Calculate precise short range Coulombic interaction energy
                   IF ( R_ij .GT. 0.0 ) THEN
                      sr_e = sr_e + e_sr(R_ij,i_x, i_y, i_z, j_x, j_y, j_z )
                   END IF
                END IF
             END IF
          END DO
       END DO
       ! Calculate layer averaged Coulombic energy from site
       disk_e = disk_e + e_disk (rho_charge(j_x), i_x, j_x)
    END DO

    ! Add short range precise energy, subtract disk averaged energies
    matrix(i_x, i_y, i_z)%esr = sr_e
    matrix(i_x, i_y, i_z)%edisk =  disk_e
    !    PRINT *, "Sr_e=",sr_e
    !    PRINT *, "Edisk=",disk_e
    !    PRINT *, "Coulomb-E=",sr_e - disk_e
    RETURN
  END SUBROUTINE coulomb_sphere

  SUBROUTINE da2l()
    ! Double-precision Array 2 Lattice
    !
    ! This subroutine copies the values of an array and copies them to the lattice
    IMPLICIT NONE
    INTEGER :: xx,yy,zz,ic

    putspin: DO ic=1,NSPIN
       ! Extract coordinates from the spin-location array
       xx = SPINLOCS(ic,1)
       yy = SPINLOCS(ic,2)
       zz = SPINLOCS(ic,3)

       ! Repopulate lattice at the appropriate locations
       socmd: IF ( (SOC_ON .EQV. .TRUE.) .OR. (MD_ON .EQV. .TRUE.) ) THEN
          ! Update lattice with external array values
          MATRIX(xx,yy,zz)%x = X(ic)
          MATRIX(xx,yy,zz)%y = Y(ic)
          MATRIX(xx,yy,zz)%z = Z(ic)
          ! Update lattice with external velocities
          MATRIX(xx,yy,zz)%vx = VX(ic)
          MATRIX(xx,yy,zz)%vy = VY(ic)
          MATRIX(xx,yy,zz)%vz = VZ(ic)
          ! Update lattice with externally calculated forces
          MATRIX(xx,yy,zz)%fx = FX(ic)
          MATRIX(xx,yy,zz)%fy = FY(ic)
          MATRIX(xx,yy,zz)%fz = FZ(ic)
       END IF socmd
       ! Update lattice with externally calculated spins
       MATRIX(xx,yy,zz)%sx = SX(ic)
       MATRIX(xx,yy,zz)%sy = SY(ic)
       MATRIX(xx,yy,zz)%sz = SZ(ic)
       socon: IF ( SOC_ON .EQV. .TRUE. ) THEN
          ! Update lattice with externally calculated SOC factors
          MATRIX(xx,yy,zz)%socx = SOCX(ic)
          MATRIX(xx,yy,zz)%socy = SOCY(ic)
          MATRIX(xx,yy,zz)%socz = SOCZ(ic)
       END IF socon
    END DO putspin
  END SUBROUTINE da2l

  subroutine DataSort (in_val)
    IMPLICIT NONE
    DOUBLE PRECISION, dimension(:,:), intent(inout):: in_val
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: val2
    integer :: in_size
    integer, ALLOCATABLE, dimension(:) :: out_index(:)
    INTEGER :: ix,ii,iix,jj,jix,nn

    in_size = SIZE(in_val,1)
    ALLOCATE( val2(in_size,SIZE(in_val,2)) )
    ALLOCATE( out_index(in_size) )
    val2 = 0

    do ix = 1, in_size
       out_index(ix) = ix
    end do

    ! Selection sort
    do ii = 1, in_size - 1
       iix = out_index(ii)
       ix = ii
       do jj = ii + 1, in_size
          jix = out_index(jj)
          ! Do your comparison here
          if (in_val(iix,1) .gt. in_val(jix,1)) then
             ! Record the smallest
             ix = jj
             iix = jix
          end if
       end do
       ! Swap
       out_index(ix) = out_index(ii)
       out_index(ii) = iix
    end do

    DO nn=1,in_size
       jj = out_index(nn)
       !      PRINT *, jj,'=',in_val(jj,:)
       val2(nn,:) = in_val(jj,:)
       !      PRINT *, nn,'=',val2(nn,:)
       !      PRINT *, "******************************************"
    END DO

    in_val = val2
    RETURN
  end subroutine DataSort

  SUBROUTINE dl2a()
    ! Double-precision Lattice 2 Array
    !
    ! This subroutine  populates an array with certain values from the lattice
    ! Values operated on are:
    ! 1) i    -> real coords
    ! 2) Vi   -> velocities
    ! 3) Fi   -> forces
    ! 4) Bi   -> HFI fields
    ! 5) Si   -> spins
    ! 6) SOCi -> spin-orbit-coupling
    IMPLICIT NONE
    INTEGER :: xx,yy,zz,ic

    ic = 0
    ! Initialize arrays
    SPINLOCS = 0

    IF ( (SOC_ON .EQV. .TRUE.) .OR. (MD_ON .EQV. .TRUE.)) THEN
       X = 0
       Y = 0
       Z = 0
       VX = 0
       VY = 0
       VZ = 0
       FX = 0
       FY = 0
       FZ = 0
    END IF

    BX = 0
    BY = 0
    BZ = 0
    SX = 0
    SY = 0
    SZ = 0

    IF ( SOC_ON .EQV. .TRUE. ) THEN
       SOCX = 0
       SOCY = 0
       SOCZ = 0
    END IF

    ! The length of the arrays should equal the number of hopping spinners in the
    ! lattice. Go through lattice and extract these values for later computation
    coord3: DO zz=1,DIMENS(3)
       coord2: DO yy=1,DIMENS(2)
          coord1: DO xx=1,DIMENS(1)
             spintest: IF ( MATRIX(xx,yy,zz)%secondary_species .NE. 0 ) THEN
                ic = ic + 1
                ! Save coordinate for later repopulation of lattice
                SPINLOCS(ic,1) = xx
                SPINLOCS(ic,2) = yy
                SPINLOCS(ic,3) = zz
                switchtest: IF ( (MD_ON .EQV. .TRUE.) .OR. (SOC_ON .EQV. .TRUE. ) ) THEN
                   ! Output real coordinates of lattice sites
                   X(ic) = MATRIX(xx,yy,zz)%x
                   Y(ic) = MATRIX(xx,yy,zz)%y
                   Z(ic) = MATRIX(xx,yy,zz)%z
                   ! Output velocities of lattice sites
                   VX(ic) = MATRIX(xx,yy,zz)%vx
                   VY(ic) = MATRIX(xx,yy,zz)%vy
                   VZ(ic) = MATRIX(xx,yy,zz)%vz
                   ! Output forces at each lattice site
                   FX(ic) = MATRIX(xx,yy,zz)%fx
                   FY(ic) = MATRIX(xx,yy,zz)%fy
                   FZ(ic) = MATRIX(xx,yy,zz)%fz
                   ! If the loop has already found all the spinners
                   ! no need to keep going. Return to calling program
                   IF ( ic .EQ. NSPIN ) RETURN
                END IF switchtest
                ! Output HFI fields at each lattice site
                BX(ic) = MATRIX(xx,yy,zz)%bx
                BY(ic) = MATRIX(xx,yy,zz)%by
                BZ(ic) = MATRIX(xx,yy,zz)%bz
                ! Output spins of each lattice site
                SX(ic) = MATRIX(xx,yy,zz)%sx
                SY(ic) = MATRIX(xx,yy,zz)%sy
                SZ(ic) = MATRIX(xx,yy,zz)%sz
                IF ( SOC_ON .EQV. .TRUE. ) THEN
                   ! Output spin-orbit-coupling at each site SOCX(ic) = MATRIX(xx,yy,zz)%socx
                   SOCY(ic) = MATRIX(xx,yy,zz)%socy
                   SOCZ(ic) = MATRIX(xx,yy,zz)%socz
                END IF
                IF ( ic .EQ. NSPIN ) RETURN
             END IF spintest
          END DO coord1
       END DO coord2
    END DO coord3
  END SUBROUTINE dl2a

  SUBROUTINE findk()
    !     COMPUTE TIME DERIVATIVES USING EQUATIONS OF MOTION (SEC. 3.3)
    !     (solve the equations of motion according to formula)
    !
    ! Note: This subroutine operates on the following data structures
    !
    ! BX  , BY  , BZ
    ! WSX , WSY , WSZ
    ! WKX , WKY , WKZ
    ! SOCX, SOCY, SOCZ
    !
    ! NB: Ask about SX vs. WSX for computation
    IMPLICIT NONE
    INTEGER :: ispin, ineib
    INTEGER :: ix,iy,iz,jx,jy,jz
    INTEGER :: ictemp
    DOUBLE PRECISION :: bexx(NSPIN),bexy(NSPIN),bexz(NSPIN)

    ! Initialize Bexi arrays
    BexX = 0.0d0
    BexY = 0.0d0
    BexZ = 0.0d0

    !     dSi/dt = Si X B
    !
    !     COMPUTE DERIVATIVES
    icloop: DO  ispin=1,NSPIN
       ! Get coordinates of cell #ispin in lattice
       CALL findcoords(ispin,ix,iy,iz)

       neibloop: do ineib=1,NNEIB
          CALL hopping(ix,iy,iz,jx,jy,jz,ineib)
          IF ( (ix .EQ. NTHICK) .AND. &
               (ineib .EQ. 5) .AND. &
               (INFINITE_HOPS .EQV. .FALSE.)) THEN
             CONTINUE
          ELSE IF ( (ix .EQ. 1) .AND. &
               (ineib .EQ. 1) .AND. &
               (INFINITE_HOPS .EQV. .FALSE.)) THEN
             CONTINUE
          ELSE
             ictemp = findspin(jx,jy,jz)
             IF ( ictemp .GT. 0 ) THEN
                BexX(ispin) = BexX(ispin) + SX(ictemp)
                BexY(ispin) = BexY(ispin) + SY(ictemp)
                BexZ(ispin) = BexZ(ispin) + SZ(ictemp)
             END IF
          END IF
       enddo neibloop

       WKX(ispin)=CHFI*(SY(ispin)*BZ(ispin)-SZ(ispin)*BY(ispin))+     &
            CJ*(SY(ispin)*BexZ(ispin)-SZ(ispin)*BexY(ispin))
       IF ( SOC_ON .EQV. .TRUE.) WKX(ispin) = WKX(ispin) - &
            CSOC*(SY(ispin)*SOCZ(ispin)-SZ(ispin)*SOCY(ispin))

       WKY(ispin)=CHFI*(SZ(ispin)*BX(ispin)-SX(ispin)*BZ(ispin))+     &
            CJ*(SZ(ispin)*BexX(ispin)-SX(ispin)*BexZ(ispin))
       IF ( SOC_ON .EQV. .TRUE. ) WKY(ispin) = WKY(ispin) - &
            CSOC*(SZ(ispin)*SOCX(ispin)-SX(ispin)*SOCZ(ispin))

       WKZ(ispin)=CHFI*(SX(ispin)*BY(ispin)-SY(ispin)*BX(ispin))+     &
            CJ*(SX(ispin)*BexY(ispin)-SY(ispin)*BexX(ispin))
       IF ( SOC_ON .EQV. .TRUE. ) WKZ(ispin) = WKZ(ispin) - &
            CSOC*(SX(ispin)*SOCY(ispin)-SY(ispin)*SOCZ(ispin))
    END DO icloop
    RETURN
  END SUBROUTINE findk

  SUBROUTINE hopping ( i_x, i_y, i_z, j_x, j_y, j_z, prob )
    !
    ! Purpose:
    !   The purpose of this  is to move a species from one site to another.
    !  It can move fr/ba/le/ri and up/down. Periodic boundary conditions are in
    !  place such that lateral motion moves to the other side of the lattice if
    !  it goes "overboard"
    !
    !  Status: 0 -> hopped
    !          1 -> could not hop
    !
    !! HOPPING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    !*****************
    ! Input and output
    !*****************
    INTEGER            , INTENT(IN)                          :: i_x, i_y, i_z
    INTEGER            , INTENT(OUT)                         :: j_x, j_y, j_z
    INTEGER            , INTENT(IN)                          :: prob

    ! Determine the direction of travel based on input number
    SELECT CASE (prob)
    CASE (1)
       j_x = i_x
       j_y = i_y
       j_z = merge(1,i_z+1,i_z.EQ.DIMENS(3))
    CASE (2)
       j_x = i_x
       j_y = i_y
       j_z = merge(DIMENS(3),i_z-1,i_z.EQ.1)
    CASE (3)
       j_x = i_x
       j_y = merge(1,i_y+1,i_y.EQ.DIMENS(2))
       j_z = i_z
    CASE (4)
       j_x = i_x
       j_y = merge(DIMENS(2),i_y-1,i_y.EQ.1)
       j_z = i_z
    CASE (5)
       j_x = merge(1,i_x + 1,((j_x .EQ. DIMENS(1)) .AND. (INFINITE_HOPS .EQV. .TRUE.)))
       j_y = i_y
       j_z = i_z
    CASE (6)
       j_x = merge(DIMENS(1),i_x - 1,((j_x .EQ. 1) .AND. (INFINITE_HOPS .EQV. .TRUE.)))
       j_y = i_y
       j_z = i_z
    END SELECT
  END SUBROUTINE hopping

  SUBROUTINE inject_electron(root,temp)
    !
    ! Purpose:
    !  The purpose of this subroutine is to inject an electron
    !  into the system at a random site on the injecting electrode,
    !  i.e., layer x=1 of the lattice.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20161004      C. Shingledecker     Original code
    !
    !! INJECT_ELECTRON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    TYPE(node), POINTER :: root, temp
    REAL :: rnd1, rnd2
    INTEGER :: i_x,i_y,i_z
    INTEGER :: j_x,j_y,j_z
    LOGICAL :: empty_site

    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Injecting electron into system!!"

    ! Initialize variables
    empty_site = .FALSE.
    i_x = 0
    i_y = 0
    i_z = 0

    ! Get ranom y and z Coordinates
    DO WHILE ( empty_site .EQV. .FALSE. )
       CALL RANDOM_NUMBER(rnd1)
       CALL RANDOM_NUMBER(rnd2)
       IF ( INITIAL_ELECS .EQV. .TRUE. ) THEN 
          i_x = 1 + FLOOR(DBLE(L+1-1)*rnd1)
          CALL RANDOM_NUMBER(rnd1)
       ELSE
          i_x = 1
       END IF
       i_y = 1 + FLOOR(DBLE(NEDGE+1-1)*rnd1)
       i_z = 1 + FLOOR(DBLE(NEDGE+1-1)*rnd2)
       IF ( MATRIX(i_x,i_y,i_z)%secondary_species .EQ. 0 ) empty_site = .TRUE.
    END DO
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "electron coordinates are:",i_x,i_y,i_z
    temp => MATRIX(i_x,i_y,i_z)
    NSPIN = NSPIN + 1

    ! Add electron notation
    temp%secondary_species = 1
    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "added secondary species"

    ! Initialize spins up
    temp%sx = 1.0d0
    temp%sy = 0.0d0
    temp%sz = 0.0d0

    ! Calculate image charge interactions
    temp%eself = e_self(i_x)

    ! Update total energy at site
    CALL sum_energies(i_x,i_y,i_z)
    CALL update_rate_and_time(i_x,i_y,i_z)

    ! Add the new electron to the tree
    CALL add_node(root,temp)
  END SUBROUTINE inject_electron

  FUNCTION e_uncorrelated (i_x)
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  E_rand_i component of the lattice site energy. The energies are
    !  uncorrelated and are sampled from a Gaussian density
    !  of states.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160920      C. Shingledecker     Original code
    !
    !! E_UNCORRELATED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION :: e_uncorrelated
    DOUBLE PRECISION :: g_E ! Gaussian random number
    DOUBLE PRECISION :: const_term
    INTEGER          :: i_x


    const_term = 1./(SIGMA_DOS*(A_NN**3)*SQRT(2*PI))
    IF ( (i_x .EQ. 1) .OR. (i_x .EQ. NTHICK) ) THEN
       e_uncorrelated = 0.0
    ELSE
       !      g_E = r8_normal_01()
       CALL RANDOM_NUMBER(g_E)
       !      PRINT *, "g_e=",g_E
       !      PRINT *, "cont_term=",const_term
       !      PRINT *, g_E/const_term
       e_uncorrelated = r8_normal_01() !(-2*(SIGMA_DOS**2)*DLOG(g_E/const_term))**0.5
       !      PRINT *, "In the function, e_un=",e_uncorrelated
    END IF
    RETURN
  END FUNCTION e_uncorrelated

  FUNCTION e_applied (i_x)
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate the
    !  E_applied_i component of the lattice site energy.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160920      C. Shingledecker     Original code
    !
    !! ENERGY_APPLIED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER :: i_x
    DOUBLE PRECISION :: e_applied

    e_applied = VC2EV*(DBLE(NTHICK - i_x)/DBLE(NTHICK - 1))*ECHARGE*VOLTAGE

  END FUNCTION e_applied

  FUNCTION e_self (i_x)
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate the
    !  E_self_i component of the lattice site energy.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160920      C. Shingledecker     Original code
    !
    !! E_SELF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION :: e_self
    DOUBLE PRECISION :: prefac
    DOUBLE PRECISION :: sum1, sum2
    INTEGER :: i_x ! Value of coord1 or lattice site
    INTEGER :: n ! Loop counter

    prefac = -1*((ECHARGE**2)/(16*PI*EPS_0*EPS_R*A_NN))
    ! Convert Volt Coulombs to eV
    prefac = prefac*VC2EV
    sum1 = 0
    sum2 = 0
    DO n=-N_IMAGES, N_IMAGES
       sum1 = sum1 + 1./(2.*i_x + 2.*DBLE(n)*DBLE(NTHICK))
       IF ( n .NE. 0 ) THEN
          sum2 = sum2 + 1./(2.*DBLE(n)*DBLE(NTHICK))
       ELSE
          sum2 = sum2 + 0.0
       END IF
    END DO
    e_self = prefac*(sum1 - sum2)
    IF ( (e_self .GT. 1.0E20) .OR. (e_self .LT. -1.0E20) ) e_self = 0.0
    RETURN
  END FUNCTION e_self

  FUNCTION e_sr (R_ij, i_x,i_y,i_z,j_x,j_y,j_z)
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  E_sr short range Coulombic lattice site energy.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160920      C. Shingledecker     Original code
    !
    !! E_SR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION :: e_sr
    INTEGER :: i_x, i_y, i_z ! Coordinates of central site
    INTEGER :: j_x, j_y, j_z ! Coordinates on lattice to calculate energy of
    INTEGER :: n
    DOUBLE PRECISION :: prefac
    DOUBLE PRECISION :: R_ij ! Distance between two points on the lattice
    DOUBLE PRECISION :: fc_ij
    DOUBLE PRECISION :: sum1, sum2

    sum1 = 0.0
    sum2 = 0.0
    prefac = ECHARGE/(4*PI*EPS_0*EPS_R)
    ! First loop to account for image charge interactions
    DO n=-N_IMAGES,N_IMAGES
       sum1 = sum1 + 1./SQRT((A_NN**2)*(DBLE(j_x + i_x + 2*n*NTHICK)**2) + &
            (A_NN**2)*(DBLE(j_y - i_y)**2) + &
            (A_NN**2)*(DBLE(j_z - i_z)**2))
       IF ( n .NE. 0 ) THEN
          sum2 = sum2 + 1./SQRT((A_NN**2)*(DBLE(j_x - i_x + 2*n*NTHICK)**2) + &
               (A_NN**2)*(DBLE(j_y - i_y)**2) + &
               (A_NN**2)*(DBLE(j_z - i_z)**2))
       ELSE
          sum2 = sum2 + 0.0
       END IF
    END DO
    ! Sum components
    fc_ij = prefac*( (1./R_ij) - sum1 + sum2 )
    e_sr = VC2EV*ECHARGE*fc_ij
  END FUNCTION e_sr

  FUNCTION e_layer (rho_charge)
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  layer averaged Coulombic energy due to carrier-carrier
    !  interactions.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160921      C. Shingledecker     Original code
    !
    !! E_LAYER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NTHICK) :: rho_charge
    DOUBLE PRECISION, DIMENSION(NTHICK) :: e_layer
    DOUBLE PRECISION, DIMENSION(NTHICK-2,NTHICK-2) :: A
    DOUBLE PRECISION, DIMENSION(NTHICK-2) :: B
    DOUBLE PRECISION                      :: a_nn_inv
    INTEGER, DIMENSION(NTHICK-2) :: IPIV
    INTEGER :: info
    INTEGER :: n

    ! Initialize values in array
    A = 0.0
    e_layer  = 0.0
    a_nn_inv = 1./A_NN

    ! Add third dimension to ρ_charge to make units work
    DO n = 1,SIZE(rho_charge,1)
       rho_charge(n) = rho_charge(n)*a_nn_inv
    END DO

    ! Make array to input
    DO n=1,NTHICK-2
       A(n,n) = 2
       IF ( n .EQ. 1 ) THEN
          A(n,2) = -1
       ELSE IF ( n .EQ. NTHICK-2 ) THEN
          A(n,n-1) = -1
       ELSE
          A(n,n-1) = -1
          A(n,n+1) = -1
       END IF
    END DO

    ! Populate the array of RHS values
    B = rho_charge(2:NTHICK-1)

    ! Solve for X using LAPACK
    ! Results will be saved in B and will have
    ! voltages that can be used to calculate layer
    ! averaged Coulombic interaction energy
    CALL sgesv(NTHICK-2,1,A,NTHICK-2,IPIV,B,NTHICK-2,info)

    ! Saver solutions to e_layer, with first and last element
    ! remaining 0
    DO n=2,NTHICK-1
       e_layer(n) = VC2EV*ECHARGE*B(n-1)
    END DO
    RETURN
  END FUNCTION e_layer

  FUNCTION e_disk (rho_charge, i_x, j_x)
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  layer averaged Coulombic energy in a spherical region
    !  considered explicitly in the E_sr, and subtract it from the
    !  layer averaged value that was previous calculated.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160921      C. Shingledecker     Original code
    !
    !! E_DISK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION :: rho_charge ! Charge density
    DOUBLE PRECISION :: e_disk ! layer averaged disk energy
    DOUBLE PRECISION :: sum1, sum2
    DOUBLE PRECISION :: sum1a, sum1b, sum2a, sum2b
    DOUBLE PRECISION :: sum1a_insides, sum2a_insides
    DOUBLE PRECISION :: prefac
    DOUBLE PRECISION :: R_ijx ! Radius of disk
    DOUBLE PRECISION :: f_disk
    INTEGER :: j_x ! current layer
    INTEGER :: i_x ! Coordinates of site at the center of the sphere
    INTEGER :: n

    prefac = rho_charge/(2.*EPS_0*EPS_R)
    sum1 = 0
    sum2 = 0
    R_ijx = SQRT((R_c**2) - (A_NN**2)*(DBLE(j_x - i_x))**2)
    DO n=-N_DISKS,N_DISKS
       sum1a_insides = (R_ijx**2) + (A_NN**2)*((DBLE(j_x + i_x + 2*n*NTHICK))**2)
       sum1a =  SQRT(sum1a_insides)
       sum1b =  A_NN*ABS(DBLE(j_x + i_x + 2*n*NTHICK))
       sum1 = sum1 + (sum1a - sum1b)
       IF ( n .NE. 0 ) THEN
          sum2a_insides = (R_ijx**2) + (A_NN**2)*((DBLE(j_x - i_x + 2*n*NTHICK)**2))
          sum2a = SQRT(sum2a_insides)
          sum2b = A_NN*ABS(DBLE(j_x - i_x + 2*n*NTHICK))
          sum2 = sum2 + (sum2a - sum2b)
       ELSE
          sum2 = sum2 + 0.0
       END IF
    END DO
    f_disk = prefac*(SQRT((R_ijx**2) + (A_NN**2)*(DBLE(j_x) - DBLE(i_x))**2) &
         - A_NN*ABS(DBLE(j_x) - DBLE(i_x)) &
         - sum1 &
         + sum2 )
    e_disk = VC2EV*ECHARGE*f_disk
    RETURN
  END FUNCTION e_disk

  SUBROUTINE sum_energies (i_x,i_y,i_z)
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate the total energy
    !  at some lattice site, i, as the sum of all of the contributing
    !  energies
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20161003      C. Shingledecker     Original code
    !
    !! SUM_ENERGIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i_x, i_y, i_z

    ASSOCIATE ( esr => MATRIX(i_x,i_y,i_z)%esr, &
         edisk    => MATRIX(i_x,i_y,i_z)%edisk, &
         eapplied => MATRIX(i_x,i_y,i_z)%eapplied, &
         erand    => MATRIX(i_x,i_y,i_z)%erand, &
         elayer   => MATRIX(i_x,i_y,i_z)%elayer, &
         eself    => MATRIX(i_x,i_y,i_z)%eself, &
         etotal   => MATRIX(i_x,i_y,i_z)%etotal)
      etotal = erand + eapplied + eself + esr + elayer - edisk
    END ASSOCIATE
  END SUBROUTINE sum_energies

  SUBROUTINE layer_current_density (current_density)
    !
    ! Purpose:
    !  The purpose of this subroutine is to calculate the linear current density
    !  of the system from one electrode to the other.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20161002      C. Shingledecker     Original code
    !
    !! LAYER_CURRENT_DENSITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(NTHICK) :: current_density
    INTEGER :: xx, yy, zz

    current_density = 0
    DO xx=1,NTHICK
       DO yy=1,NEDGE
          DO zz=1,NEDGE
             current_density(xx) = current_density(xx) + MATRIX(xx,yy,zz)%J_i
          END DO
       END DO
    END DO
    OPEN(UNIT=LAYERUNIT,FILE='current_density.wsv',POSITION="APPEND")
    WRITE(LAYERUNIT,*) current_density
    CLOSE(LAYERUNIT)
  END SUBROUTINE layer_current_density

  FUNCTION miller_abrahams ( delta_e_ij )
    !
    ! Purpose:
    !  The purpose of this function is to calculate the
    !  Miller-Abrahams hopping rate between site i and
    !  site j
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160922      C. Shingledecker     Original code
    !
    !! MILLER_ABRAHAMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    DOUBLE PRECISION :: miller_abrahams ! The Miller-Abrahams hopping rate
    !    DOUBLE PRECISION :: R_ij ! Distance between sites i and j
    DOUBLE PRECISION :: delta_e_ij ! Energy difference between sites i and j
    CHARACTER(len=1) :: cont

    IF ( DEBUG .EQV. .TRUE. ) THEN
       PRINT *, "ν0 =",NU_0 ! Trial frequency
       PRINT *, "α =",ALPHA !
       PRINT *, "a_nn=",A_NN ! Fixed distance to nearest neighbor
       PRINT *, "K_b=",K_B ! Boltzmann constant
       PRINT *, "T_Kelvin=",T_KELVIN ! Temperature in Kelvin
       PRINT *, "Δeij=",delta_e_ij ! Energy difference between lattice sites
    END IF

    IF ( delta_e_ij  .GT. 0.0 ) THEN
       miller_abrahams = NU_0*EXP((-2*ALPHA*A_NN))*EXP((-1*delta_e_ij)/(K_B*T_KELVIN))
    ELSE
       miller_abrahams = NU_0*EXP(-2*ALPHA*A_NN)
    END IF
    RETURN
  END FUNCTION miller_abrahams

  SUBROUTINE new_lattice(temp)
    ! This subroutine reproduces the functionality of lattice.f, but is compatible with
    ! the KMC code.
    !
    ! 20161022 -> Original Code
    !
    ! Christopher N. Shingledecker
    IMPLICIT NONE
    INTEGER :: xx,yy,zz ! Coordinates
    INTEGER :: ic ! Unique identifier of each lattice site
    TYPE(node), POINTER :: temp ! temporary holder for matrix node

    !******************************************************************************
    ! Create the matrix
    !
    ! NOTE: Incorporates features of lattice.f
    !
    !******************************************************************************
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, 'In main, the dimens are:',dimens
    ALLOCATE ( MATRIX( DIMENS(1),DIMENS(2),DIMENS(3) ) )

    ! Initialize the matrix
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "In main, matrix has been initialized to 0, now assigning -1 to O2"
    ic = 0
    coord3: DO zz=1,dimens(3)
       coord2: DO yy=1,dimens(2)
          coord1: DO xx=1,dimens(1)
             ic = ic + 1
             temp => MATRIX(xx,yy,zz)
             CALL init_node(temp)
             !**********************************
             ! Populate values from spin/md code
             !**********************************
             ! Assign unique identifier
             temp%ic = ic
             ! Integer coordinates in lattice
             temp%ix = xx
             temp%iy = yy
             temp%iz = zz
             ! Real coordinates in lattice
             temp%x = DBLE(xx)*DISTNN
             temp%y = DBLE(yy)*DISTNN
             temp%z = DBLE(zz)*DISTNN
             ! Random initial velocities
             temp%vx = r8_normal_01()*MBSCALE
             temp%vy = r8_normal_01()*MBSCALE
             temp%vz = r8_normal_01()*MBSCALE
             IF ( NO_HOP .EQV. .TRUE. ) THEN
                ! Initial spin states (all pointed up: i.e. x=1)
                temp%secondary_species = 1
                temp%sx = 1.0d0
             ELSE
                temp%sx = 0.0d0
             END IF
             temp%sy = 0.0d0
             temp%sz = 0.0d0
             ! Initial HFI fields
             temp%bx = r8_normal_01()
             temp%by = r8_normal_01()
             temp%bz = r8_normal_01()
             ! Populate values from electron transport/kmc code
             IF ( (xx .NE. 1) .AND. (xx .NE. dimens(1)) ) THEN
                temp%primary_species = 1 ! Here we assign -1 to sites occupied by the major bulk species
             ELSE
                temp%primary_species = 0 ! Electrode value
             END IF
             temp%erand = e_uncorrelated(xx)
             temp%eapplied = e_applied(xx)
             CALL sum_energies(xx,yy,zz)
             IF ( (ASSOCIATED(temp%before)) .OR. (ASSOCIATED(temp%after)) ) THEN
                PRINT *, "Uh oh, something is still associated"
                CALL EXIT()
             END IF
          END DO coord1
       END DO coord2
    END DO coord3
    ! Save the number of elements in the lattice
    NELEMS = ic
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Matrix has been fully initialized"
  END SUBROUTINE new_lattice

  SUBROUTINE PREDCOR()
    !     FOURTH-ORDER PREDICTOR-CORRECTOR METHOD (SEC. 3.3)
    !     Adams-Bashforth predictor Adams-Moulton Corrector
    !     SEE NUMERICAL ANALYSIS by Burden & Faires pg. 279
    IMPLICIT NONE
    integer*4 ispin
    real*8 C5,C9,PC9,C19,C37,C59,C55
    DATA C5,C9,PC9,C19,C37,C59,C55&
         /-5.d0,-9.d0,9.d0,19.d0,37.d0,-59.d0,55.d0/

    !     PREDICT
    DO ispin=1,nspin
       WSX(ispin)=SX(ispin)+DT24*                                     &
            (C9*SX0(ispin)+C37*SX1(ispin)+C59*SX2(ispin)              &
            +C55*SX3(ispin))
       WSY(ispin)=SY(ispin)+DT24*                                     &
            (C9*SY0(ispin)+C37*SY1(ispin)+C59*SY2(ispin)              &
            +C55*SY3(ispin))
       WSZ(ispin)=SZ(ispin)+DT24*                                     &
            (C9*SZ0(ispin)+C37*SZ1(ispin)+C59*SZ2(ispin)              &
            +C55*SZ3(ispin))
    ENDDO

    CALL FINDK()

    DO ispin=1,nspin
       SX(ispin)=SX(ispin)+DT24*                                      &
            (SX1(ispin)+C5*SX2(ispin)+C19*SX3(ispin)+PC9*WKX(ispin))
       SY(ispin)=SY(ispin)+DT24*                                      &
            (SY1(ispin)+C5*SY2(ispin)+C19*SY3(ispin)+PC9*WKY(ispin))
       SZ(ispin)=SZ(ispin)+DT24*                                      &
            (SZ1(ispin)+C5*SZ2(ispin)+C19*SZ3(ispin)+PC9*WKZ(ispin))
       SX0(ispin)=SX1(ispin)
       SY0(ispin)=SY1(ispin)
       SZ0(ispin)=SZ1(ispin)
       SX1(ispin)=SX2(ispin)
       SY1(ispin)=SY2(ispin)
       SZ1(ispin)=SZ2(ispin)
       SX2(ispin)=SX3(ispin)
       SY2(ispin)=SY3(ispin)
       SZ2(ispin)=SZ3(ispin)
       WSX(ispin)=SX(ispin)
       WSY(ispin)=SY(ispin)
       WSZ(ispin)=SZ(ispin)
    ENDDO

    CALL FINDK()

    DO ispin=1,nspin
       SX3(ispin)=WKX(ispin)
       SY3(ispin)=WKY(ispin)
       SZ3(ispin)=WKZ(ispin)
    ENDDO

    RETURN
  END SUBROUTINE PREDCOR

  SUBROUTINE randomize_hfi()
    IMPLICIT NONE
    INTEGER :: xx, yy, zz
    DO zz=1,DIMENS(3)
       DO yy=1,DIMENS(2)
          DO xx=1,DIMENS(1)
             MATRIX(xx,yy,zz)%bx = r8_normal_01()
             MATRIX(xx,yy,zz)%by = r8_normal_01()
             MATRIX(xx,yy,zz)%bz = r8_normal_01()
          END DO
       END DO
    END DO
  END SUBROUTINE randomize_hfi

  FUNCTION rc_int_dist (i,j,k) !
    ! Purpose:
    !  The purpose of this subroutine is to find the number of
    !  integer hops one must traverse to exceed the critical
    !  radius, R_c, in the matrix.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160922      C. Shingledecker     Original code
    !
    !! RC_INT_DIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER, DIMENSION(6) :: rc_int_dist
    INTEGER :: i, j, k ! Coordinates on lattice to calculate energy of
    DOUBLE PRECISION :: R_ij ! Distance between two points on the lattice

    rc_int_dist = 0
    ! Hop outwards in one direction to find dimensions of square to analyze.
    ! Note: Here we approximate the square in later x as being defined by
    ! one value, l/2, where l is the length of one side of the square.
    ! In principal, though, one could also use xmin, xmax, jmin, jmax, etc.

    ! Find z_max
    DO
       IF ( k+rc_int_dist(1) .GE. dimens(3) ) EXIT
       rc_int_dist(1) = rc_int_dist(1) + 1
       R_ij = ABS(matrix(i,j,k)%z - matrix(i,j,k+rc_int_dist(1))%z)
       IF ( R_ij .GE. R_c ) EXIT
    END DO

    ! Find z_min
    DO
       IF ( k+rc_int_dist(2) .LE. 1 ) EXIT
       rc_int_dist(2) = rc_int_dist(2) - 1
       R_ij = ABS(matrix(i,j,k)%z - matrix(i,j,k+rc_int_dist(2))%z)
       IF ( R_ij .GE. R_c ) EXIT
    END DO

    ! Find y_max
    DO
       IF ( j+rc_int_dist(3) .GE. dimens(2) ) EXIT
       rc_int_dist(3) = rc_int_dist(3) + 1
       R_ij = ABS(matrix(i,j,k)%y - matrix(i,j+rc_int_dist(3),k)%y)
       IF ( R_ij .GE. R_c ) EXIT
    END DO

    ! Find y_min
    DO
       IF ( j+rc_int_dist(4) .LE. 1 ) EXIT
       rc_int_dist(4) = rc_int_dist(4) - 1
       R_ij = ABS(matrix(i,j,k)%y - matrix(i,j+rc_int_dist(4),k)%y)
       IF ( R_ij .GE. R_c ) EXIT
    END DO

    ! Find x_max
    DO
       IF ( i+rc_int_dist(5) .GE. dimens(1) ) EXIT
       rc_int_dist(5) = rc_int_dist(5) + 1
       R_ij = ABS(matrix(i,j,k)%x - matrix(i+rc_int_dist(5),j,k)%x)
       IF ( R_ij .GE. R_c ) EXIT
    END DO

    ! Find x_min
    DO
       IF ( i+rc_int_dist(6) .LE. 1 ) EXIT
       rc_int_dist(6) = rc_int_dist(6) - 1
       R_ij = ABS(matrix(i,j,k)%x - matrix(i+rc_int_dist(6),j,k)%x)
       IF ( R_ij .GE. R_c ) EXIT
    END DO
    RETURN
  END FUNCTION rc_int_dist

  SUBROUTINE rungek4()
    !     FOURTH-ORDER RUNGE-KUTTA METHOD (SEC. 3.3)
    !     SEE NUMERICAL ANALYSIS by Burden & Faires pg. 259
    !
    ! Note: This subroutine operates on the following global data structures
    !
    ! SX , SY , SZ
    ! SX0, SY0, SZ0
    ! WSX, WSY, WSZ
    ! WKX, WKY, WKZ
    !
    IMPLICIT NONE
    INTEGER :: ispin
    !     FIND K1

    DO 100 ispin=1,NSPIN
       WSX(ispin)=SX0(ispin)
       WSY(ispin)=SY0(ispin)
       WSZ(ispin)=SZ0(ispin)
100 END DO

    CALL FINDK()

    DO 200 ispin=1,NSPIN
       SX(ispin)=SX0(ispin)+DT6*WKX(ispin)
       SY(ispin)=SY0(ispin)+DT6*WKY(ispin)
       SZ(ispin)=SZ0(ispin)+DT6*WKZ(ispin)
200 END DO

    !     FIND K2

    DO 300 ispin=1,NSPIN
       WSX(ispin)=SX0(ispin)+DT2*WKX(ispin)
       WSY(ispin)=SY0(ispin)+DT2*WKY(ispin)
       WSZ(ispin)=SZ0(ispin)+DT2*WKZ(ispin)
300 END DO

    CALL FINDK()

    DO 400 ispin=1,NSPIN
       SX(ispin)=SX(ispin)+DT3*WKX(ispin)
       SY(ispin)=SY(ispin)+DT3*WKY(ispin)
       SZ(ispin)=SZ(ispin)+DT3*WKZ(ispin)
400 END DO

    !      FIND K3

    DO 500 ispin=1,NSPIN
       WSX(ispin)=SX0(ispin)+DT2*WKX(ispin)
       WSY(ispin)=SY0(ispin)+DT2*WKY(ispin)
       WSZ(ispin)=SZ0(ispin)+DT2*WKZ(ispin)
500 END DO

    CALL FINDK()

    DO 600 ispin=1,NSPIN
       SX(ispin)=SX(ispin)+DT3*WKX(ispin)
       SY(ispin)=SY(ispin)+DT3*WKY(ispin)
       SZ(ispin)=SZ(ispin)+DT3*WKZ(ispin)
600 END DO

    !      FIND K4

    DO 700 ispin=1,NSPIN
       WSX(ispin)=SX0(ispin)+DT*WKX(ispin)
       WSY(ispin)=SY0(ispin)+DT*WKY(ispin)
       WSZ(ispin)=SZ0(ispin)+DT*WKZ(ispin)
700 END DO

    CALL FINDK()

    DO 800 ispin=1,NSPIN
       SX(ispin)=SX(ispin)+DT6*WKX(ispin)
       SY(ispin)=SY(ispin)+DT6*WKY(ispin)
       SZ(ispin)=SZ(ispin)+DT6*WKZ(ispin)
800 END DO

    RETURN
  END SUBROUTINE rungek4

  SUBROUTINE update_lattice (root, temp, prevNode, nextNode, hopcount)
    !
    ! Purpose:
    !  The purpose of this subroutine is to assign a value of
    !  E_interact, i.e. a Coulombic energy value, to each site
    !  in the lattice. This subroutine is called after each
    !  electron hop, i.e. when the relative position of charge
    !  carriers changes in the model and a new value of the
    !  energy needs to be computed.
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160925      C. Shingledecker     Original code
    !
    !! UPDATE_LATTICE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER :: error
    INTEGER :: hopcount
    DOUBLE PRECISION, DIMENSION(NTHICK) :: rho_charge
    DOUBLE PRECISION, DIMENSION(NTHICK) :: layer_energy
    TYPE (node), POINTER :: root
    TYPE (node), POINTER :: temp
    TYPE (node), POINTER :: prevNode, nextNode
    INTEGER :: xx, yy, zz

    NULLIFY( temp, prevNode, nextNode )

    ! Step 1: Count the charge carriers in each layer of the lattice
    rho_charge = charge_density()

    ! Step 2: Calculate the layer-averaged Coulombic interaction energy
    !         and assign this value as E_interact in each lattice layer,
    !         if the new value is different from the old.
    layer_energy = e_layer(rho_charge)

    ! Step 3: If R_c is not 0, go to each site in the lattice and calculate the
    !         explicit short range charge interaction in a sphere of radius R_c
    !         around each site. Add this new value, E_sr, to the E_interact
    !         value and subtract the E_disk value corresponding to the layer value
    !         that would otherwise be double-counted.
    DO xx=1,NTHICK
       DO yy=1,NEDGE
          DO zz=1,NEDGE
             matrix(xx,yy,zz)%elayer = layer_energy(xx)
             ! IF ( DEBUG .EQV. .TRUE. ) PRINT *, "In update lattice at:",xx,yy,zz,"fixed e=",&
             !                                    matrix(xx,yy,zz)%e_fixed
             ! IF ( DEBUG .EQV. .TRUE. ) PRINT *, "Now calling coulomb sphere"
             IF ( R_c .GT. 0 ) CALL coulomb_sphere(xx,yy,zz,rho_charge)
             ! Update total energy at lattice site
             CALL sum_energies(xx,yy,zz)
             ! If lattice site contains an electron, update hopping rates
             IF ( MATRIX(xx,yy,zz)%secondary_species .NE. 0 ) THEN
                temp => MATRIX(xx,yy,zz)
                ! Delete node from tree, since wait time will be updated
                !            PRINT *, "Calling delete node"
                CALL delete_node(root, temp, prevNode, nextNode, error)
                temp => MATRIX(xx,yy,zz)
                ! Update KMC values
                ! PRINT *, "Calling update rate and time"
                CALL update_rate_and_time(xx, yy, zz)
                ! Re-add to tree with new time
                ! PRINT *, "Calling update add_node"
                CALL add_node(root, temp)
             END IF
             ! Randomize HFI fields
             matrix(xx,yy,zz)%bx = r8_normal_01()
             matrix(xx,yy,zz)%by = r8_normal_01()
             matrix(xx,yy,zz)%bz = r8_normal_01()
             !***End of Inner Loop***
          END DO
       END DO
    END DO
  END SUBROUTINE update_lattice

  SUBROUTINE update_rate_and_time (i_x, i_y, i_z)
    !
    ! Purpose:
    !  The purpose of this subroutine is to update the hopping rate
    !  and time for a particular site based on the ΔE values around
    !  it. Based on these, one is selected from it's relative magnatice
    !  and a new hopping time is calculated
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20160925      C. Shingledecker     Original code
    !
    !! UPDATE_RATE_AND_TIME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE parameters
    IMPLICIT NONE
    INTEGER :: i_x, i_y, i_z
    INTEGER :: j_x, j_y, j_z
    INTEGER :: n
    INTEGER :: temp_ind
    DOUBLE PRECISION, DIMENSION(6,3) :: rates
    DOUBLE PRECISION :: rand
    DOUBLE PRECISION :: tau
    DOUBLE PRECISION :: temp_rate
    DOUBLE PRECISION :: delta_e_ij
    DOUBLE PRECISION :: rate

    ! Initialize rates to 0
    rates = 0.0
    temp_rate = 0.0

    ! Calculate rate to nearest neighbors
    ! The cases should correspond to the cases
    ! in the hopping subroutine
    CALL sum_energies(i_x,i_y,i_z)
    DO n=1,6
       rates(n,2) = n
       CALL hopping ( i_x, i_y, i_z, j_x, j_y, j_z, n )
       !      PRINT *, "j_coords=",j_x,j_y,j_z
       IF ( (j_x .NE. 0) .AND. (j_y .NE. 0) .AND. (j_z .NE. 0)) THEN
          CALL sum_energies(j_x,j_y,j_z)
          delta_e_ij = MATRIX(j_x,j_y,j_z)%etotal - MATRIX(i_x,i_y,i_z)%etotal
          IF ( MATRIX(j_x,j_y,j_z)%secondary_species.EQ.0 ) THEN
             rates(n,1) = miller_abrahams(delta_e_ij)
          ELSE
             rates(n,1) = 0.0
          END IF
       ELSE
          rates(n,1) = 0.0
       END IF
       rates(n,3) = rates(n,1)
    END DO

    IF ( i_x .NE. 1 ) THEN
       CALL DataSort(rates)
       CALL normalize_probability_and_select_index(rates,temp_ind)
       ! PRINT *, "Just after normalize_probability_and_select_index, temp_ind=",temp_ind

       ! Call a new random number to get waiting time, τ, for hop
       rate = 0
       DO n=1,6
          IF ( INT(rates(n,2)) .EQ. temp_ind ) rate = rates(n,3)
       END DO
    ELSE
       rate = rates(5,3)
       temp_ind = 5
    END IF

    ! PRINT *, "The rate is:",rate
    CALL RANDOM_NUMBER(rand)
    IF ( rate .EQ. 0.0 ) THEN
       DO n=1,6
          PRINT *, rates(n,:)
          IF ( rates(n,3) .GT. rate ) THEN
             rate = rates(n,3)
             temp_ind = n
          END IF
       END DO
       IF ( rate .EQ. 0.0 ) THEN
          MATRIX(i_x,i_y,i_z)%hop_dir = 0
          MATRIX(i_x,i_y,i_z)%tau =  TIME + TAU
          MATRIX(i_x,i_y,i_z)%omega = 0
          RETURN
       END IF
    END IF

    taucalc: IF ( FIXEDSTEP .EQV. .TRUE. ) THEN
       tau = FIXEDTAU
    ELSE
       tau = -1*(LOG(rand)/rate)
    END IF taucalc

    IF ( DEBUG .EQV. .TRUE. ) PRINT *, "τ=",tau

    ! Update matrix with new rate, hopping direction, and waiting time (τ)
    MATRIX(i_x,i_y,i_z)%tau = tau + TIME
    MATRIX(i_x,i_y,i_z)%omega = rate
    MATRIX(i_x,i_y,i_z)%hop_dir = temp_ind
  END SUBROUTINE update_rate_and_time

  SUBROUTINE write_matrix_dataframe(i_x,i_y,i_z)
    !
    ! Purpose:
    !  The purpose of this subroutine is to go through the matrix
    !  and write the values of each lattice site
    !
    ! Documentation:
    !  DATE          PROGRAMMER           DESCRIPTION
    !  ========      ==========           ===========
    !  20161003      C. Shingledecker     Original code
    !
    !! WRITE_MATRIX_DATAFRAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE
    INTEGER :: i_x, i_y, i_z
    CHARACTER(len=80) :: fmt

    OPEN(UNIT=MATRIXUNIT,FILE='matrix_printout.csv',POSITION="APPEND")
    fmt = "(I2,A1,I2,A1,I2,A1,ES15.6,A1,A9)"
    IF ( matrix(i_x,i_y,i_z)%secondary_species .NE. 0 ) THEN
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%tau,',',"tau"
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%omega,',',"omega"
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%eself,',',"E_self"
    END IF

    WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%erand,',',"E_rand"
    WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%eapplied,',',"E_applied"

    IF ( matrix(i_x,i_y,i_z)%elayer .NE. 0 ) THEN
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%elayer,',',"E_layer"
    END IF

    IF ( matrix(i_x,i_y,i_z)%esr .NE. 0 ) THEN
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%esr,',',"E_sr"
    END IF

    IF ( matrix(i_x,i_y,i_z)%edisk .NE. 0 ) THEN
       WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%edisk,',',"E_disk"
    END IF

    WRITE(MATRIXUNIT,fmt) i_x,',',i_y,',',i_z,',',matrix(i_x,i_y,i_z)%etotal,',',"E_total"
    CLOSE(MATRIXUNIT)
  END SUBROUTINE write_matrix_dataframe

  SUBROUTINE findcoords(ic,ix,iy,iz)
    ! This subroutine finds the coordinates of a site in a lattice, given a unique cell
    ! identifier
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ic
    INTEGER, INTENT(OUT) :: ix, iy, iz

    ix = 0
    iy = 0
    iz = 0

    ! 1. Calculate z
    iz = CEILING(REAL(ic)/(REAL(NTHICK)*REAL(NEDGE)))
    ! 2. Calculate y
    iy = CEILING((REAL(ic)-(iz-1)*(REAL(NTHICK)*REAL(NEDGE)))/REAL(NTHICK))
    ! 3. Calculate x
    ix = (ic - (iz-1)*(REAL(NTHICK)*REAL(NEDGE))) - (iy-1)*NTHICK
    RETURN
  END SUBROUTINE findcoords

  SUBROUTINE normalize_probability_and_select_index(array,temp_ind)
    IMPLICIT NONE
    DOUBLE PRECISION :: array(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: temp_array(:,:)
    DOUBLE PRECISION :: prev_val
    DOUBLE PRECISION :: rn
    DOUBLE PRECISION :: rate_sum
    INTEGER :: nn, n
    INTEGER, INTENT(OUT) :: temp_ind
    INTEGER :: not_zero, temp_counter

    prev_val = 0.0
    not_zero = 0
    temp_ind = 0
    prev_val = 0.0

    rate_sum = SUM(array(:,1))
    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) THEN
       DO n=1,6
          !      rates(n,1) = rates(n,1)/rate_sum
          PRINT *, array(n,:)
       END DO
    END IF

    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) PRINT *, "After division by rates sum:"
    array(:,1) = array(:,1)/rate_sum
    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) THEN
       DO n=1,6
          PRINT *, array(n,:)
       END DO
    END IF

    ! Note: This only works on a list that has been sorted in ascending order
    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) PRINT *, "After normalization:"
    DO nn=2,SIZE(array,1)+1
       array(nn-1,1) = array(nn-1,1) + prev_val
       prev_val = array(nn-1,1)
       IF ( array(nn-1,1) .GT. 0.0 ) not_zero = not_zero + 1
       IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) PRINT *, array(nn-1,:)
    END DO

    IF ( not_zero .EQ. 0 ) THEN
       PRINT *, "Not_zero = 0"
       CALL EXIT()
    END IF

    ALLOCATE(temp_array(not_zero,SIZE(array,2)))
    temp_array = 0
    temp_counter = 1
    DO n=1,SIZE(array,1)
       IF ( array(n,1) .GT. 0.0 ) THEN
          temp_array(temp_counter,:) = array(n,:)
          !        PRINT *, "temp_array(",temp_counter,",:)=",temp_array(temp_counter,:)
          temp_counter = temp_counter + 1
       END IF
    END DO

    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) THEN
       PRINT *, "The not_zero array is:"
       DO n=1,not_zero
          PRINT *, temp_array(n,:)
       END DO
    END IF

    CALL RANDOM_NUMBER(rn)
    temp_ind = INT(temp_array(not_zero,2))
    DO nn=1,SIZE(temp_array,1)
       IF ( ( rn .GT. prev_val ) .AND. ( rn .LE. temp_array(nn,1) ) ) THEN
          temp_ind = INT(temp_array(nn,2))
       END IF
       prev_val = temp_array(nn,1)
    END DO

    IF ( (DEBUG .EQV. .TRUE.) .AND. (VERBOSE .EQV. .TRUE.) ) PRINT *, "The temp_ind=",temp_ind
    ! IF ( temp_ind .EQ. 0 ) THEN
    !   DO n=1,SIZE(array,1)
    !     IF ( (array(n,2) .EQ. 5.0) .AND. (array(n,1) .NE. 0.0) ) THEN
    !       temp_ind = n
    !     ELSE IF ( array(n,1) .NE. 0.0 ) THEN
    !       temp_ind = n
    !     END IF
    !   END DO
    ! END IF
  END SUBROUTINE normalize_probability_and_select_index

  SUBROUTINE transfer(ix,iy,iz,jx,jy,jz,root,temp,prevNode,nextNode)
    INTEGER   , INTENT(IN)          :: ix, iy, iz ! Initial coords of hopper
    INTEGER   , INTENT(IN)          :: jx, jy, jz ! Final coords of hopper
    INTEGER                         :: error
    TYPE(node)            , POINTER :: root, temp, prevNode, nextNode

    error = 0

    !***************************************************************************
    ! Now that the electron has moved, update energies in lattice for next hop
    !***************************************************************************
    collectOrHop: IF ( (jx .NE. NTHICK) .OR. (INFINITE_HOPS .EQV. .TRUE.) ) THEN
       ! Repopulate new lattice site with electron informatino
       temp                   => MATRIX(jx,jy,jz)
       temp%secondary_species = 1
       temp%eself             = e_self(jx)
       temp%sx                = MATRIX(ix,iy,iz)%sx
       temp%sy                = MATRIX(ix,iy,iz)%sy
       temp%sz                = MATRIX(ix,iy,iz)%sz
       temp%wsx               = MATRIX(ix,iy,iz)%wsx
       temp%wsy               = MATRIX(ix,iy,iz)%wsy
       temp%wsz               = MATRIX(ix,iy,iz)%wsz
       temp%wkx               = MATRIX(ix,iy,iz)%wkx
       temp%wky               = MATRIX(ix,iy,iz)%wky
       temp%wkz               = MATRIX(ix,iy,iz)%wkz
       CALL sum_energies(jx,jy,jz)
       CALL update_rate_and_time(jx,jy,jz)
       CALL add_node(root,temp)
    ELSE
       PRINT *, "Electron arriving at collecting electrode"
       NSPIN = NSPIN - 1
    END IF collectOrHop

    !***************************************************************************
    ! Delete from tree and wipe old lattice site (electron no longer there)
    !***************************************************************************
    temp => MATRIX(ix,iy,iz)
    CALL delete_node(root,temp,prevNode,nextNode,error)
    temp => MATRIX(ix,iy,iz)
    CALL wipe_node(temp)
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Node deleted"
  END SUBROUTINE transfer

  FUNCTION findspin(ix,iy,iz)
    INTEGER             :: findspin
    INTEGER, INTENT(IN) :: ix,iy,iz
    INTEGER             :: coords(3)
    INTEGER             :: ii

    findspin  = 0

    coords(1) = ix
    coords(2) = iy
    coords(3) = iz

    ! See if input coords are in spinlocs array
    DO ii = 1,SIZE(SPINLOCS,1)
       IF ( ALL(ABS(SPINLOCS(1,:) - coords) .EQ. 0) ) findspin = ii
    END DO

    RETURN
  END FUNCTION findspin
END MODULE subroutines
