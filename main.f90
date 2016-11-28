PROGRAM main
  USE spincalc
  USE subroutines
  USE parameters
  USE typedefs
  USE bsimple
  IMPLICIT NONE
  !******************************************************************************
  ! Main/KMC Variables
  !******************************************************************************
  INTEGER                   :: xx,yy,zz       ! Counters
  INTEGER                   :: count_num      ! Number of abundance file
  INTEGER                   :: time_check     ! DEBUGGING VAR
  INTEGER                   :: error
  INTEGER                   :: i_x, i_y, i_z
  INTEGER                   :: j_x, j_y, j_z
  INTEGER                   :: hop_count
  DOUBLE PRECISION          :: time_diff      !DEBUGGING VAR
  DOUBLE PRECISION          :: t1,t2
  DOUBLE PRECISION          :: cpu_total
  DOUBLE PRECISION          :: cpu_max_time
  DOUBLE PRECISION          :: J_ix(NTHICK)
  TYPE (node)     , POINTER :: root, temp
  TYPE (node)     , POINTER :: prevNode, nextNode
  INTEGER                   :: itime, ispin, ineib, ihfi
  INTEGER                   :: ii

  !To enable debugging outputs, set debug to true
  IF ( VERBOSE .EQV. .TRUE. ) THEN
     PRINT *, "*************************"
     PRINT *, "***STARTING SIMULATION***"
     PRINT *, "*************************"
  END IF

  ! Replace old output files
  OPEN(UNIT=1, FILE="pxt.wsv", STATUS="REPLACE")
  OPEN(UNIT=2, FILE="pxt_deriv.wsv", STATUS="REPLACE")
  OPEN(UNIT=3, FILE="new_kthop.wsv", STATUS="REPLACE")
  OPEN(UNIT=4, FILE="new_kt.wsv", STATUS="REPLACE")
  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
  CLOSE(4)

  time_check = 0
  time_diff  = 0
  cpu_total  = 0
  t2         = 0
  t1         = 0

  ! Allocate pointers
  ALLOCATE ( root, temp, prevNode, nextNode )

  ! Nullify pointers
  NULLIFY ( root, temp, prevNode, nextNode )

  ! Initialize time step between abundance checks
  t1 = 0

  ! Initialize count_num
  count_num = 1
  cpu_max_time = 100.0
  cpu_total = 0
  TIME = 0.0

  ! Build the lattice
  CALL new_lattice(temp)

  ! Populate the lattice with electrons, if appropriate
  IF ( (INITIAL_ELECS .EQV. .TRUE.) .AND. (NO_HOP .EQV. .FALSE. ) ) THEN
     DO ii = 1,EMAX
        initelec: IF ( NSPIN .LT. EMAX ) THEN
           IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Now adding electron"
           CALL inject_electron(root,temp)
           IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Finished adding electron"
        END IF initelec
     END DO
  END IF


  hop_count = 0
  IF ( NO_HOP .EQV. .TRUE. ) THEN
     NEXTHOP = MAXTIME
     NSPIN   = L*L*L
  END IF

  !******************************************************************************
  ! Begin the simulation
  !******************************************************************************
  PRINT *, "Now starting main loop"
  mainloop: DO WHILE ( (TIME+DT) .LT. MAXTIME )
     !IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Loop #",hop_count, "NSPIN=",NSPIN,"Time=",TIME

     nohop1: IF ( NO_HOP .EQV. .FALSE. ) THEN
        !***************************************************************************
        ! Find min time event
        !***************************************************************************
        findmin: IF ( ASSOCIATED(root) ) THEN
           CALL find_min(root, temp)
        ELSE
           CALL inject_electron(root,temp)
        END IF findmin

        !***************************************************************************
        ! Save next hop time
        !***************************************************************************
        NEXTHOP = temp%tau
        IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "NEXTHOP=",NEXTHOP - TIME
     END IF nohop1

     !***************************************************************************
     ! Relax Spins
     !***************************************************************************
     relaxspin: IF ((TIME+(3*DT)) .LT. NEXTHOP )  THEN
        IF ( NSPIN .GT. 0 ) THEN
           IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "NSPIN=",NSPIN,"CALLING spinrelax"
           CALL spinrelax()
        END IF
     END IF relaxspin

     nohop2: IF ( NO_HOP .EQV. .FALSE. ) THEN
        !***************************************************************************
        ! Update time and coordinates for hop
        !***************************************************************************
        TIME = NEXTHOP
        i_x = temp%ix
        i_y = temp%iy
        i_z = temp%iz
        j_x = 0
        j_y = 0
        j_z = 0

        !***************************************************************************
        ! Update site current density values
        !***************************************************************************
        IF ( MATRIX(i_x,i_y,i_z)%hop_dir .EQ. 5 ) THEN           ! Increment f_i at site
           MATRIX(i_x,i_y,i_z)%f_i = MATRIX(i_x,i_y,i_z)%f_i + 1 ! Update current density (C/s)
        ELSE IF ( MATRIX(i_x,i_y,i_z)%hop_dir .EQ. 6 ) THEN      ! Increment b_i at site
           MATRIX(i_x,i_y,i_z)%b_i = matrix(i_x,i_y,i_z)%b_i + 1 ! Update current density (C/s)
        ELSE IF ( MATRIX(i_x,i_y,i_z)%hop_dir .EQ. 0 ) THEN      ! If hopdir = 0, do nothing
           GOTO 900
        END IF

        ! Calculate new site current density
        MATRIX(i_x,i_y,i_z)%J_i = (ECHARGE*DBLE(MATRIX(i_x,i_y,i_z)%f_i - &
             MATRIX(i_x,i_y,i_z)%b_i))/TIME

        IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "J_i=",MATRIX(i_x,i_y,i_z)%J_i

        !***************************************************************************
        ! Move electron
        !***************************************************************************
        hop_count = hop_count + 1
        CALL hopping(i_x,i_y,i_z,j_x,j_y,j_z,temp%hop_dir)
        IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "j_coords are", j_x,j_y,j_z
        IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Hopdir=",temp%hop_dir,&
             ':',i_x,i_y,i_z, &
             "->",j_x,j_y,j_z

        !***************************************************************************
        ! Move the electron, i.e. transfer unique vals. to new lattice site and
        ! wipe old lattice site of same vals.
        !***************************************************************************
        IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Now deleting node"
        CALL transfer(i_x,i_y,i_z,j_x,j_y,j_z,root,temp,prevNode,nextNode)

        !***************************************************************************
        ! Update lattice energies with new charge positions
        !***************************************************************************
        IF ( FIXEDSTEP .EQV. .TRUE. ) THEN
           CALL update_rate_and_time(j_x,j_y,j_z)
        ELSE
           CALL update_lattice(root, temp, prevNode, nextNode, hop_count)
           CALL layer_current_density(J_ix)
        END IF


        !***************************************************************************
        ! If hop_dir = 0, Go here, i.e. do not hop or update energies
        !***************************************************************************
900     CONTINUE

        !***************************************************************************
        ! Add electron to system, if appropriate
        !***************************************************************************
        IF ( INITIAL_ELECS .EQV. .FALSE. ) THEN
           addelec: IF ( NSPIN .LT. EMAX ) THEN
              IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Now adding electron"
              CALL inject_electron(root,temp)
              IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Finished adding electron"
           END IF addelec
        END IF
     END IF nohop2

     !***************************************************************************
     ! Analytics: Time, etc.
     !***************************************************************************
     time_check = time_check + 1
     CALL CPU_TIME(t2)
     cpu_total = cpu_total + (t2-t1)
     time_diff = TIME
     t1 = t2
     !     IF ( VERBOSE .EQV. .TRUE. ) PRINT *, &
     !          "Time=",TIME,&
     !          "CPU total=",cpu_total,&
     !          "time diff=",time_diff
  END DO mainloop


  IF ( VERBOSE .EQV. .TRUE. ) THEN
     PRINT *, "*****************"
     PRINT *, "ENDING SIMULATION"
     PRINT *, "*****************"
  END IF
END PROGRAM main
