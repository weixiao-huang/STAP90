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
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
! .                                                                   . 
! .   P R O G R A M                                                   . 
! .        TO CALCULATE ISOPARAMETRIC HEXAHEDRAL ELEMENT STIFFNESS    . 
! .        MATRIX                                                     . 
! .                                                                   . 
! .  - - INPUT VARIABLES - -                                          . 
! .        NEL       = NUMBER OF ELEMENT                              . 
! .        NINT      = GAUSS NUMERICAL INTEGRATION ORDER              . 
! .        E          = YOUNG'S MODULUS                               . 
! .        PR        = POISSON'S RATIO                                .
! .        S(24,24)    = STORAGE FOR STIFFNESS MATRIX                 . 
! .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    . 
! .                                                                   . 
! .  - - OUTPUT - -                                                   . 
! .        S(24,24)    = CALCULATED STIFFNESS MATRIX                  . 
! .                                                                   . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

SUBROUTINE Hexahedral
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the Hexahedral element subroutine    .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 2*NUMMAT*ITWO +25*NUME + 24*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: PR(NUMMAT)
! N103: LM(24,NUME)
! N104: XYZ(24,NUME)
! N105: MTAP(NUME)

  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+24*NUME
  N105=N104+24*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL exahedral (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE Hexahedral


SUBROUTINE exahedral (ID,X,Y,Z,U,MHT,E,PR,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   Hexahedral element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(3,NUMNP),LM(24,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),PR(NPAR(3)),  &
             XYZ(24,NPAR(2)),U(NEQ),UE(24)                           
  REAL(8) :: S(24,24),D(6,6)
  
  INTEGER :: NPAR1, NUME, NUMMAT, ND,P1,P2,P3,P4,P5,P6,P7,P8, L, N, I,J,K
  INTEGER :: MTYPE, IPRINT                                       
  REAL(8) :: STR(6), P(6)
  REAL(8) :: GP(2), WGT(2),B(6,24),detJ
  REAL(8),parameter:: pi=3.1415926535D0
  GP =[-0.5773502692 , 0.5773502692]                                            !两点高斯积分
  WGT=[   1          ,   1         ]
  
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 
  ND = 24

! Read and generate element information
  IF (IND .EQ. 1) THEN

        WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                    ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                    '     EQ.1, TRUSS ELEMENTS',/,      &
                    '     EQ.2, 4Q ELEMENTS',/,  &
                    '     EQ.3, 3T ELEMENTS',/,      &
                    '     EQ.4, 8H ELEMENTS',/,      &
                    '     EQ.5, BEAM ELEMENTS',/,  &
                    '     EQ.6, PLATE ELEMENTS',/,      &
                    '     EQ.7, SHELL ELEMENTS',//,      &
                    ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                 ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                 ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                 4 (' .'),'( NPAR(3) ) . . =',I5,//)")  NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S      POISSION ',/,  &
                   ' NUMBER     MODULUS         RATIO  ',/,  &
                   '               E             PR    ')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,F10.0,f10.5)') N,E(N),PR(N)                      ! Read material information  
        WRITE (IOUT,"(I5,4X,E12.5,2X,E12.5)") N,E(N),PR(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE    NODE    NODE    NODE    NODE    NODE    NODE    NODE    MATERIAL',/,   &
                      ' NUMBER-N     P1      P2      P3      P4      P5      P6      P7      P8     SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(10I5)') N,P1,P2,P3,P4,P5,P6,P7,P8,MTYPE  ! Read in element information

!       Save element information
        XYZ(1,N)=X(P1)     ! Coordinates of the element's first node
        XYZ(2,N)=Y(P1)
        XYZ(3,N)=Z(P1)
        XYZ(4,N)=X(P2)     ! Coordinates of the element's second node
        XYZ(5,N)=Y(P2)      
        XYZ(6,N)=Z(P2)
        XYZ(7,N)=X(P3)     ! Coordinates of the element's third node
        XYZ(8,N)=Y(P3)
        XYZ(9,N)=Z(P3)
        XYZ(10,N)=X(P4)     ! Coordinates of the element's fourth node
        XYZ(11,N)=Y(P4)
        XYZ(12,N)=Z(P4)
        XYZ(13,N)=X(P5)     ! Coordinates of the element's fifth node
        XYZ(14,N)=Y(P5)
        XYZ(15,N)=Z(P5)
        XYZ(16,N)=X(P6)     ! Coordinates of the element's sixth node
        XYZ(17,N)=Y(P6)
        XYZ(18,N)=Z(P6)
        XYZ(19,N)=X(P7)     ! Coordinates of the element's seventh node
        XYZ(20,N)=Y(P7)
        XYZ(21,N)=Z(P7)
        XYZ(22,N)=X(P8)     ! Coordinates of the element's eighth node
        XYZ(23,N)=Y(P8)
        XYZ(24,N)=Z(P8)
        MATP(N)=MTYPE       ! Material type
        
        DO L=1,24
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)   =ID(L,P1)     ! Connectivity matrix
           LM(L+3,N) =ID(L,P2) 
           LM(L+6,N) =ID(L,P3) 
           LM(L+9,N) =ID(L,P4)
           LM(L+12,N)=ID(L,P5)
           LM(L+15,N)=ID(L,P6)
           LM(L+18,N)=ID(L,P7)
           LM(L+21,N)=ID(L,P8)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,7X,I5)") N,P1,P2,P3,P4,P5,P6,P7,P8,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN
    
      DO N=1,NUME
        MTYPE=MATP(N)
        
        D=0
        D(1,1)=1-PR(MTYPE)
        D(1,2)=PR(MTYPE)
        D(1,3)=PR(MTYPE)
        D(2,1)=PR(MTYPE)
        D(2,2)=1-PR(MTYPE)
        D(2,3)=PR(MTYPE)
        D(3,1)=PR(MTYPE)
        D(3,2)=PR(MTYPE)
        D(3,3)=1-PR(MTYPE)
        D(4,4)=0.5-PR(MTYPE)
        D(5,5)=0.5-PR(MTYPE)
        D(6,6)=0.5-PR(MTYPE)
        D=D*(E(MTYPE)/((1+PR(MTYPE))*(1-2*PR(MTYPE))))
  
        S=0
   
        DO I=1,2                                                      
            DO J=1,2
                DO K=1,2
                    CAll Bmatr(GP(I),GP(J),GP(K),XYZ(1,N),B,detJ)

                    S=S+WGT(I)*WGT(J)*WGT(K)*matmul(matmul(transpose(B),D),B)*detJ
                END DO
            END DO
        END DO
 
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
                                           '  ELEMENT',4X,'GAUSS POINT',15X,'StressXX',9X,'StressYY',9X,'StressZZ',9X,'StressXY',9X,'StressYZ',9X,'StressXZ',/,&
                                          '  NUMBER')") NG
        MTYPE=MATP(N)
        do L=0,7    
            I=LM(3*L+1,N)
            IF (I.GT.0) then
                UE(3*L+1)=U(I)
            else
                UE(3*L+1)=0
            end if       
            
            J=LM(3*L+2,N)            
            if (J.GT.0) then
                UE(3*L+2) = U(J)
            else
                UE(3*L+2) = 0
            end if
            
            K=LM(3*L+3,N)            
            if (K.GT.0) then
                UE(3*L+3) = U(K)
            else
                UE(3*L+3) = 0
            end if

        end do
        
        DO I=1,2                                                      
            DO J=1,2
                DO K=1,2
                    CALL Bmatr(GP(I),GP(J),GP(K),XYZ(1,N),B,detJ)
                    STR = matmul(B,UE)
                    P   = matmul(D,STR)
                    write (IOUT,"(I5,5X,f6.3,2X,f6.3,2X,f6.3,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6)")N,GP(I),GP(J),GP(K),P(1),P(2),P(3),P(4),P(5),P(6)
                END DO
            end do
        end do

     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE exahedral

    
    
subroutine Bmatr(xi,eta,zeta,XY,B,detJ)
    implicit none
    real(8) :: xi,eta,zeta,XY(24),loca(8,3),B(6,24),detJ,GNxez(3,8),GNxyz(3,8),J(3,3),JINV(3,3),DUM
    integer:: K,I

    GNxez(1,1)=-(1-eta)*(1+zeta)
    GNxez(1,2)=(1-eta)*(1+zeta)
    GNxez(1,3)=(1-eta)*(1-zeta)
    GNxez(1,4)=-(1-eta)*(1-zeta)
    GNxez(1,5)=-(1+eta)*(1+zeta)
    GNxez(1,6)=(1+eta)*(1+zeta)
    GNxez(1,7)=(1+eta)*(1-zeta)
    GNxez(1,8)=-(1+eta)*(1-zeta)
    GNxez(2,1)=-(1-xi)*(1+zeta)
    GNxez(2,2)=-(1+xi)*(1+zeta)
    GNxez(2,3)=-(1+xi)*(1-zeta)
    GNxez(2,4)=-(1-xi)*(1-zeta)
    GNxez(2,5)=(1-xi)*(1+zeta)
    GNxez(2,6)=(1+xi)*(1+zeta)
    GNxez(2,7)=(1+xi)*(1-zeta)
    GNxez(2,8)=(1-xi)*(1-zeta)
    GNxez(3,1)=(1-xi)*(1-eta)
    GNxez(3,2)=(1+xi)*(1-eta)
    GNxez(3,3)=-(1+xi)*(1-eta)
    GNxez(3,4)=-(1-xi)*(1-eta)
    GNxez(3,5)=(1-xi)*(1+eta)
    GNxez(3,6)=(1+xi)*(1+eta)
    GNxez(3,7)=-(1+xi)*(1+eta)
    GNxez(3,8)=-(1-xi)*(1+eta)
    GNxez=0.125*GNxez

    DO K=0,7
        loca(K+1,1)=XY(3*K+1)
        loca(K+1,2)=XY(3*K+2)
        loca(K+1,3)=XY(3*K+3)
    END DO
        
    J=matmul(GNxez,loca)
    detJ=J(1,1)*J(2,2)*J(3,3)+J(1,2)*J(2,3)*J(3,1)+J(1,3)*J(2,1)*J(3,2)-J(1,3)*J(2,2)*J(3,1)-J(1,2)*J(2,1)*J(3,3)-J(1,1)*J(2,3)*J(3,2)
    DUM=1./detJ
    JINV(1,1)=DUM*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
    JINV(1,2)=DUM*(J(1,3)*J(3,2)-J(1,2)*J(3,3))
    JINV(1,3)=DUM*(J(1,2)*J(2,3)-J(1,3)*J(2,2))
    JINV(2,1)=DUM*(J(2,3)*J(3,1)-J(2,1)*J(3,3))
    JINV(2,2)=DUM*(J(1,1)*J(3,3)-J(1,3)*J(3,1))
    JINV(2,3)=DUM*(J(1,3)*J(2,1)-J(1,1)*J(2,3))
    JINV(3,1)=DUM*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
    JINV(3,2)=DUM*(J(1,2)*J(3,1)-J(1,1)*J(3,2))
    JINV(3,3)=DUM*(J(1,1)*J(2,2)-J(1,2)*J(2,1))

    GNxyz=matmul(JINV,GNxez)
    
    DO I=1,6
        DO K=1,24
            B(I,K)=0
        END DO
    END DO
    DO K=0,7
        B(1,3*K+1)=GNxyz(1,K+1)
        B(2,3*K+2)=GNxyz(2,K+1)
        B(3,3*K+3)=GNxyz(3,K+1)
        B(4,3*K+2)=GNxyz(3,K+1)
        B(4,3*K+3)=GNxyz(2,K+1)
        B(5,3*K+1)=GNxyz(3,K+1)
        B(5,3*K+3)=GNxyz(1,K+1)
        B(6,3*K+1)=GNxyz(2,K+1)
        B(6,3*K+2)=GNxyz(1,K+1)
    END DO
        
end subroutine Bmatr