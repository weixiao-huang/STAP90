! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                            S T A P 9 0                                .
! .                                                                       .
! .     AN IN-CORE SOLUTION STATIC ANALYSIS PROGRAM IN FORTRAN 90         .
! .     Adapted from STAP (KJ Bath, FORTRAN IV) for teaching purpose      .
! .                                                                       .
! .     Xiong Zhang, (2013)                                               .
! .     Computational Dynamics Group, School of Aerospace                 .
! .     Tsinghua Univerity                                                .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

SUBROUTINE ELCAL
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To loop over all element groups for reading,                    .
! .   generating and storing the element data                         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: N, I

  REWIND IELMNT
  WRITE (IOUT,"(//,' E L E M E N T   G R O U P   D A T A',//)")

! Loop over all element groups

  DO N=1,NUMEG
     IF (N.NE.1) WRITE (IOUT,'(1X)')

     READ (IIN,'(10I5)') NPAR

     CALL ELEMNT

     IF (MIDEST.GT.MAXEST) MAXEST=MIDEST

     WRITE (IELMNT) MIDEST,NPAR,(A(I),I=NFIRST,NLAST)

  END DO

  RETURN

END SUBROUTINE ELCAL


SUBROUTINE ELEMNT
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call the appropriate element subroutine                      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS

  IMPLICIT NONE
  INTEGER :: NPAR1

  NPAR1=NPAR(1)

  IF (NPAR1 == 1) THEN
     CALL TRUSS
  ELSE IF (NPAR1 == 2) THEN    ! Quadrilateral Elements
     CALL QuadrElem
  ELSE IF (NPAR1 == 3) THEN    ! Triangle Elements
     CALL TRIELMT
  ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
     CALL Hexahedral
  ELSE IF (NPAR1 == 5) THEN    ! Beam Elements
     CALL BEAM
  ELSE IF (NPAR1 == 6) THEN    ! Plate Elements
     CALL PLATE
  ELSE IF (NPAR1 == 7) THEN    ! Shell Elements
     CALL SHELL
  ELSE
!    Other element types would be called here, identifying each
!    element type by a different NPAR(1) parameter
  END IF

  RETURN
END SUBROUTINE ELEMNT


SUBROUTINE STRESS (AA)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To call the element subroutine for the calculation of stresses  .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS, ONLY : IELMNT, NG, NUMEST, NPAR, NUMEG

  IMPLICIT NONE
  REAL :: AA(*)
  INTEGER N, I

! Loop over all element groups

  REWIND IELMNT

  DO N=1,NUMEG
     NG=N

     READ (IELMNT) NUMEST,NPAR,(AA(I),I=1,NUMEST)

     CALL ELEMNT
  END DO

  RETURN
END subroutine STRESS

    
subroutine EIGEN
  USE GLOBALS
  USE MEMALLOCATE
  IMPLICIT NONE
    
    INTEGER :: NROOT, NC, NNC
    INTEGER :: NITEM,IFSS,IFPR,NSTIF, NWM
    
    REAL(8) :: TT(NEQ), W(NEQ), RTOL
    REAL(8),ALLOCATABLE :: RV(:,:), EIGV(:), AR(:), BR(:)
    REAL(8),ALLOCATABLE :: VEC(:,:), D(:), RTOLV(:), BUP(:), BLO(:)
    REAL(8),ALLOCATABLE :: BUPC(:), Q(:,:)
    
    DATA NITEM/16/, IFSS/1/, IFPR/1/, NSTIF/16/
    DATA RTOL/1.0E-6/
    
    NROOT = NPAR(6)
    NC = MIN(2*NROOT, NROOT+8, NEQ)
    NNC = NC*(NC+1)/2
    NWM = NWK

    IF (NROOT > NEQ) STOP "*** ERROR *** eigenvalue is lager than NEQ."
    
    ALLOCATE (RV(NEQ,NC),EIGV(NC),AR(NNC),BR(NNC))
    ALLOCATE (VEC(NC,NC),D(NC),RTOLV(NC),BUP(NC),BLO(NC))
    ALLOCATE (BUPC(NC),Q(NEQ,NC))
    
    IF (NPAR(5) == 1) THEN          ! 子空间迭代法
        CALL SSPACE90(DA(NP(3)),DA(NP(13)),IA(NP(2)),RV,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ, &
            NEQ+1,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)
    ELSE IF (NPAR(5) == 2) THEN     ! Lanczos法
        CALL LANCZOS(DA(NP(3)),DA(NP(13)),IA(NP(2)),RV,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ, &
            NEQ+1,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT,Q)
    ELSE
        WRITE(*,"(/,'NOT CALCULATE THE EIGENVALUES AND EIGENVECTORS.')") 
    END IF
  
  return
end subroutine EIGEN