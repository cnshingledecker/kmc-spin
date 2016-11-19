MODULE spincalc
  USE subroutines
  USE parameters
  USE typedefs
  USE bsimple
CONTAINS
  SUBROUTINE spinrelax()
    IMPLICIT NONE
    !******************************************************************************
    ! Spin/MD Variables
    !******************************************************************************
    INTEGER                   :: itime, ispin, ineib, ihfi
    DOUBLE PRECISION          :: t0,t1,t2,t3
    DOUBLE PRECISION          :: t
    DOUBLE PRECISION          :: sxave, syave, szave
    DOUBLE PRECISION          :: wkxave
    DOUBLE PRECISION          :: kt, kthop

    !To enable debugging outputs, set debug to true
    IF ( VERBOSE .EQV. .TRUE. ) THEN
       PRINT *, "********************************"
       PRINT *, "***STARTING SPIN CALCULATIONS***"
       PRINT *, "********************************"
    END IF

    OPEN(UNIT=1, FILE="pxt.wsv"      , STATUS="UNKNOWN",POSITION="APPEND")
    OPEN(UNIT=2, FILE="pxt_deriv.wsv", STATUS="UNKNOWN",POSITION="APPEND")
    OPEN(UNIT=3, FILE="new_kt.wsv"   , STATUS="UNKNOWN",POSITION="APPEND")
    OPEN(UNIT=4, FILE="new_kthop.wsv", STATUS="UNKNOWN",POSITION="APPEND")
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "Opened files"

    t3 = 0
    t2 = 0
    t1 = 0
    t0 = TIME

    !******************************************************************************
    ! SiO input in Runge-Kutta 4 subroutine
    !******************************************************************************
    ! 0) Extract spin values from lattice and save them to Si
    !******************************************************************************
    ALLOCATE(SPINLOCS(NSPIN,3))

    IF ( (SOC_ON .EQV. .TRUE.) .OR. (MD_ON .EQV. .TRUE.) ) THEN
       ALLOCATE(X(NSPIN),Y(NSPIN),Z(NSPIN))
       ALLOCATE(VX(NSPIN),VY(NSPIN),VZ(NSPIN))
       ALLOCATE(FX(NSPIN),FY(NSPIN),FZ(NSPIN))
    END IF

    ALLOCATE(BX(NSPIN),BY(NSPIN),BZ(NSPIN))
    ALLOCATE(SX(NSPIN),SY(NSPIN),SZ(NSPIN))
    ALLOCATE(SX0(NSPIN),SY0(NSPIN),SZ0(NSPIN))
    ALLOCATE(SX1(NSPIN),SY1(NSPIN),SZ1(NSPIN))
    ALLOCATE(SX2(NSPIN),SY2(NSPIN),SZ2(NSPIN))
    ALLOCATE(SX3(NSPIN),SY3(NSPIN),SZ3(NSPIN))
    ALLOCATE(WSX(NSPIN),WSY(NSPIN),WSZ(NSPIN))
    ALLOCATE(WKX(NSPIN),WKY(NSPIN),WKZ(NSPIN))

    IF ( SOC_ON .EQV. .TRUE. ) THEN
       ALLOCATE(SOCX(NSPIN),SOCY(NSPIN),SOCZ(NSPIN))
    END IF

    CALL dl2a()

    !******************************************************************************
    ! 1) Put existing spin values in Si0 and Si3
    !******************************************************************************
    DO 100 ispin=1,NSPIN
       SX0(ispin) = SX(ispin)
       SY0(ispin) = SY(ispin)
       SZ0(ispin) = SZ(ispin)
       SX3(ispin) = SX0(ispin)
       SY3(ispin) = SY0(ispin)
       SZ3(ispin) = SZ0(ispin)
100 END DO

    !******************************************************************************
    ! 2) Call Runge-Kutta algorithm and update time
    !******************************************************************************
    CALL rungek4()
    TIME = TIME + DT
    t1 = TIME

    !******************************************************************************
    ! 3) Send new Si value to Si1 and also to Si0
    !******************************************************************************
    DO 200 ispin=1,NSPIN
       SX1(ispin) = SX(ispin)
       SY1(ispin) = SY(ispin)
       SZ1(ispin) = SZ(ispin)
       SX0(ispin) = SX(ispin)
       SY0(ispin) = SY(ispin)
       SZ0(ispin) = SZ(ispin)
200 END DO

    !******************************************************************************
    ! 4) Call Runge-Kutta algorithm and update time
    !******************************************************************************
    CALL rungek4()
    TIME = TIME + DT
    t2 = TIME

    !******************************************************************************
    ! 5) Send new Si value to Si2 and also to Si0
    !******************************************************************************
    DO 300 ispin=1,NSPIN
       SX2(ispin) = SX(ispin)
       SY2(ispin) = SY(ispin)
       SZ2(ispin) = SZ(ispin)
       SX0(ispin) = SX(ispin)
       SY0(ispin) = SY(ispin)
       SZ0(ispin) = SZ(ispin)
300 END DO

    !******************************************************************************
    ! 6) Call Runge-Kutta algorithm and update time
    !******************************************************************************
    CALL rungek4()
    TIME = TIME + DT
    t3 = TIME

    !******************************************************************************
    ! 7) Si3 (that holds original values) -> Si0
    !    Si -> Si3
    !    Si0 -> WSi (to calculate derivatives)
    !******************************************************************************
    DO 400 ispin=1,NSPIN
       SX0(ispin) = SX3(ispin)
       SY0(ispin) = SY3(ispin)
       SZ0(ispin) = SZ3(ispin)
       SX3(ispin) = SX(ispin)
       SY3(ispin) = SY(ispin)
       SZ3(ispin) = SZ(ispin)
       WSX(ispin) = SX0(ispin)
       WSY(ispin) = SY0(ispin)
       WSZ(ispin) = SZ0(ispin)
400 END DO


    kt = 1./3. + 2./3.*(1-(t0*dlt)**2.0)*EXP(-0.5*t0*t0*dlt*dlt)
    kthop = DEXP(-2.d0*(dlt/khop)**2.0*(DEXP(-khop*t0)-1.d0+khop*t0))
    sxave = SUM(SX0)/DBLE(NSPIN)
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "sxave0=",sxave
    WRITE(1,*) t0, sxave
    WRITE(3,*) t0, kt
    WRITE(4,*) t0, kthop

    kt = 1./3. + 2./3.*(1-(t1*dlt)**2.0)*EXP(-0.5*t1*t1*dlt*dlt)
    kthop = DEXP(-2.d0*(dlt/khop)**2.0*(DEXP(-khop*t1)-1.d0+khop*t1))
    sxave = SUM(SX1)/DBLE(NSPIN)
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "sxave1=",sxave
    WRITE(1,*) t1, sxave
    WRITE(3,*) t1, kt
    WRITE(4,*) t1, kthop

    kt = 1./3. + 2./3.*(1-(t2*dlt)**2.0)*EXP(-0.5*t2*t2*dlt*dlt)
    kthop = DEXP(-2.d0*(dlt/khop)**2.0*(DEXP(-khop*t2)-1.d0+khop*t2))
    sxave = SUM(SX2)/DBLE(NSPIN)
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "sxave2=",sxave
    WRITE(1,*) t2, sxave
    WRITE(3,*) t2, kt
    WRITE(4,*) t2, kthop

    sxave = SUM(SX3)/DBLE(NSPIN)
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "sxave3=",sxave
    WRITE(1,*) t3, sxave

    ! Now we have Si0, Si1, Si2, Si3
    ! Compute time-derivatives using equations of motion
    CALL findk()

    ! Compute derivative of Si0 and put it in its place, then perform the same
    ! for Si1,Si2, and Si3.
    ! The original value of Si3 is kept in Si (which is the current value this is Si3
    !                                          at the end of the algorithm)
    DO 500 ispin=1,NSPIN
       SX0(ispin) = WKX(ispin)
       SY0(ispin) = WKY(ispin)
       SZ0(ispin) = WKZ(ispin)
       WSX(ispin) = SX1(ispin)
       WSY(ispin) = SY1(ispin)
       WSZ(ispin) = SZ1(ispin)
500 END DO

    CALL findk()

    DO 600 ispin=1,NSPIN
       SX1(ispin) = WKX(ispin)
       SY1(ispin) = WKY(ispin)
       SZ1(ispin) = WKZ(ispin)
       WSX(ispin) = SX2(ispin)
       WSY(ispin) = SY2(ispin)
       WSZ(ispin) = SZ2(ispin)
600 END DO

    CALL findk()

    DO 700 ispin=1,NSPIN
       SX2(ispin) = WKX(ispin)
       SY2(ispin) = WKY(ispin)
       SZ2(ispin) = WKZ(ispin)
       WSX(ispin) = SX3(ispin)
       WSY(ispin) = SY3(ispin)
       WSZ(ispin) = SZ3(ispin)
700 END DO

    CALL findk()

    DO 800 ispin=1,NSPIN
       SX3(ispin) = WKX(ispin)
       SY3(ispin) = WKY(ispin)
       SZ3(ispin) = WKZ(ispin)
800 END DO

    ! Calculate derivative of 1st time-step
    wkxave = SUM(SX0)/DBLE(NSPIN)
    WRITE(2,*) t0, wkxave
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "wkxave0=",wkxave

    ! Calculate derivative of 2nd time-step
    wkxave = SUM(SX1)/DBLE(NSPIN)
    WRITE(2,*) t1, wkxave
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "wkxave1=",wkxave

    ! Calculate derivative of 3rd time-step
    wkxave = SUM(SX2)/DBLE(NSPIN)
    WRITE(2,*) t2, wkxave
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "wkxave2=",wkxave

    ! Calculate derivative of 4th time-step
    wkxave = SUM(SX3)/DBLE(NSPIN)
    WRITE(2,*) t3, wkxave
    IF ( VERBOSE .EQV. .TRUE. ) PRINT *, "wkxave3=",wkxave

    !******************************************************************************
    ! Main Loop
    !******************************************************************************
    spinloop: DO WHILE (TIME .LT. NEXTHOP)
       TIME = TIME + DT

       CALL predcor()

       sxave = 0.d0
       DO ispin=1,NSPIN
          sxave = sxave + SX(ispin)
       END DO
       sxave = sxave/DBLE(NSPIN)

       ! Calculate analytical results
       kt = 1./3. + 2./3.*(1-(TIME*dlt)**2.0)*EXP(-0.5*TIME*TIME*dlt*dlt)
       kthop = DEXP(-2.d0*(dlt/khop)**2.0*(DEXP(-khop*TIME)-1.d0+khop*TIME))

       WRITE(1,*) TIME, sxave
       WRITE(3,*) TIME, kt
       WRITE(4,*) TIME, kthop
       PRINT *, "Time=",TIME, "Sxave=",sxave, "kt=",kt,"kthop=",kthop

       wkxave = 0.d0
       DO ispin=1,NSPIN
          wkxave = wkxave + WKX(ispin)
       END DO
       wkxave = wkxave/DBLE(NSPIN)

       WRITE(2,*) TIME, wkxave
       IF ( TIME .GE. MAXTIME ) EXIT
    END DO spinloop

    ! Close files
900 CONTINUE
    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)

    !***************************************************************************
    ! Repopulate lattice with new spins and other values
    !***************************************************************************
    CALL da2l()

    DEALLOCATE(SPINLOCS)

    IF ( (SOC_ON .EQV. .TRUE.) .OR. (MD_ON .EQV. .TRUE.) ) THEN
       DEALLOCATE(VX,VY,VZ)
       DEALLOCATE(FX,FY,FZ)
    END IF

    DEALLOCATE(BX,BY,BZ)
    DEALLOCATE(SX,SY,SZ)
    DEALLOCATE(SX0,SY0,SZ0)
    DEALLOCATE(SX1,SY1,SZ1)
    DEALLOCATE(SX2,SY2,SZ2)
    DEALLOCATE(SX3,SY3,SZ3)
    DEALLOCATE(WSX,WSY,WSZ)
    DEALLOCATE(WKX,WKY,WKZ)

    IF ( SOC_ON .EQV. .TRUE. ) THEN
       DEALLOCATE(SOCX,SOCY,SOCZ)
    END IF

    IF ( VERBOSE .EQV. .TRUE. ) THEN
       PRINT *, "***********************"
       PRINT *, "ENDING SPIN CALCULATION"
       PRINT *, "***********************"
    END IF
  END SUBROUTINE spinrelax
END MODULE spincalc
