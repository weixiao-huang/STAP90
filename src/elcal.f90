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
!     CALL PLATE
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
    
    INTEGER, PARAMETER :: NC=2
    INTEGER, PARAMETER :: NNC=NC*(NC+1)/2
    
    INTEGER :: NROOT,NITEM,IFSS,IFPR,NSTIF, NWM

    REAL(8) :: RV(NEQ,NC),EIGV(NC),TT(NEQ),W(NEQ),AR(NNC),BR(NNC)
    REAL(8) :: VEC(NC,NC),D(NC),RTOLV(NC),BUP(NC),BLO(NC)
    REAL(8) :: BUPC(NC),RTOL
    REAL(8) :: Q(NEQ,NC)

    DATA NROOT/2/, NITEM/16/, IFSS/1/, IFPR/1/, NSTIF/16/
    DATA RTOL/1.0E-6/

    NWM = NWK
    
    IF (NPAR(5) == 0) THEN
        CALL SSPACE90(DA(NP(3)),DA(NP(13)),IA(NP(2)),RV,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ, &
            NEQ+1,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT)
    ELSE IF (NPAR(5) == 1) THEN
        CALL LANCZOS(DA(NP(3)),DA(NP(13)),IA(NP(2)),RV,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,BUPC,NEQ, &
            NEQ+1,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT,Q)
    ELSE
        STOP "*** ERROR *** Invalid NPAR(5) value."
    END IF
  
  return
end subroutine EIGEN