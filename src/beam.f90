SUBROUTINE BEAM
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the beam element subroutine          .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106,N107

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 4*NUMMAT*ITWO + 7*NUME + 6*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: AREA(NUMMAT)
! N103: LM(6,NUME)
! N104: XYZ(6,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+6*NUME
  N105=N104+6*NUME*ITWO
  N106=N105+NUME
  N107=N106+NUMMAT*ITWO
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL EAM (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N106),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE BEAM


SUBROUTINE EAM (ID,X,Y,Z,U,MHT,E,AREA,INERTIAL,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   beam element subroutine                                         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(3,NUMNP),LM(6,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),AREA(NPAR(3)),  &
            INERTIAL(NPAR(3)), XYZ(6,NPAR(2)),U(NEQ),T(NUMNP)                            ! TÎª½Ç¶Ètheta
  REAL(8) :: S(6,6),ST(6),D(3)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I, J, L, N
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: XL2, XL, SQRT,XM, XX, YY, STR, P
  REAL(8) :: EAA, EII 
  
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=6

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, QUADR ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.5, BEAM ELEMENTS ',//,     &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND PROPETIES ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL    DENSITY',/,  &
                   ' NUMBER     MODULUS',10X,'AREA         VALUE',/,  &
                   15 X,'E',14X,'A',15X,'I')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,3F10.0)') N,E(N),AREA(N),INERTIAL(N)                      ! INERTIAL  ¹ßÐÔ¾Ø
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,2X,E12.5)") N,E(N),AREA(N),INERTIAL(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(5I5)') N,I,J,MTYPE  ! Read in element information

!       Save element information
        XYZ(1,N)=X(I)  ! Coordinates of the element's left node
        XYZ(2,N)=Y(I)
        XYZ(3,N)=Z(I)
        
        XYZ(4,N)=X(J)  ! Coordinates of the element's right node
        XYZ(5,N)=Y(J)
        XYZ(6,N)=Z(J)
        
        MATP(N)=MTYPE  ! Material type

        DO L=1,6
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)=ID(L,I)     ! Connectivity matrix
           LM(L+3,N)=ID(L,J)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,7X,I5)") N,I,J,MTYPE
        WRITE (10,"(I5,4X,I5)") I,J
        
     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)

        XL2=0.
        DO L=1,3
           D(L)=XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO
        XL=SQRT(XL2)   ! Length of element N
        
        EAA=E(MTYPE)*AREA(MTYPE)    !EA
        EII=INERTIAL(MTYPE)*E(MTYPE)      !IA
        
        DO I=1,6
            DO J=1,6
                S(I,J)=0
            END DO
        END DO
        
        S(1,1)=EAA/XL
        S(1,4)=-EAA/XL
        S(4,1)=-EAA/XL
        S(4,4)=EAA/XL

        S(2,2)=12*EII/XL/XL/XL
        S(2,5)=-12*EII/XL/XL/XL        
        S(5,2)=-12*EII/XL/XL/XL 
        S(5,5)=12*EII/XL/XL/XL 
        
        S(2,3)=6*EII/XL/XL
        S(2,6)=6*EII/XL/XL
        S(5,3)=-6*EII/XL/XL
        S(5,6)=-6*EII/XL/XL
        
        S(3,2)=6*EII/XL/XL
        S(3,5)=-6*EII/XL/XL
        S(6,2)=6*EII/XL/XL
        S(6,5)=-6*EII/XL/XL
        
        S(3,3)=4*EII/XL
        S(3,6)=2*EII/XL
        S(6,3)=2*EII/XL
        S(6,6)=4*EII/XL
              
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',13X,'FORCE',12X,'STRESS',/,'  NUMBER')") NG
        MTYPE=MATP(N)

        P=0
        STR=0
        
        XL2=0.
        DO L=1,3
           D(L)=XYZ(L,N) - XYZ(L+3,N)
           XL2=XL2 + D(L)*D(L)
        END DO
        XL=SQRT(XL2)   ! Length of element N
        
        EAA=E(MTYPE)*AREA(MTYPE)    !EA
        EII=INERTIAL(MTYPE)*E(MTYPE)      !IA
        
        WRITE (IOUT,"(1X,I5,11X,E13.6,4X,E13.6)") N,P,STR
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE EAM
