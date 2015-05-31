! (NOTE: NOT HAVE THE DENSITY FOR THE MASS MATRIX. NOW DENSITY IS EQUAL TO 1e3)
SUBROUTINE TRIELMT6
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the triangle element subroutine      .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106,N107

  NUME = NPAR(2)
  NUMMAT = NPAR(3)
! ITYPE = NPAR(4)
! EIGEN = NPAR(5)     ! advanced elements

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 3*NUMMAT*ITWO +7*NUME + 6*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: PR(NUMMAT)
! N103: LM(6,NUME)
! N104: XY(6,NUME)
! N105: MTAP(NUME)
! N106: THIC(NUMMAT)
  N101 = NFIRST
  N102 = N101 + NUMMAT * ITWO
  N103 = N102 + NUMMAT * ITWO
  N104 = N103 + 6 * NUME
  N105 = N104 + 6 * NUME * ITWO
  N106 = N105 + NUME
  N107 = N106 + NUMMAT * ITWO
  NLAST = N107

  MIDEST = NLAST - NFIRST

  CALL triangle6 (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),IA(NP(5)),   &
                 A(N101),A(N102),A(N106),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE TRIELMT6


SUBROUTINE triangle6(ID,X,Y,U,MHT,E,PR,THIC,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   triangle element subroutine                                     .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    USE GLOBALS
    USE MEMALLOCATE

    IMPLICIT NONE
    INTEGER :: ID(3,NUMNP), LM(6,NPAR(2)), MATP(NPAR(2)), MHT(NEQ)
    REAL(8) :: X(NUMNP), Y(NUMNP), E(NPAR(3)), &
               PR(NPAR(3)), THIC(NPAR(3)), &
               XY(6,NPAR(2)), U(NEQ)
    INTEGER :: NPAR1, NUME, NUMMAT, ND, II(3)
    INTEGER :: I, J, L, N
    INTEGER :: MTYPE, IPRINT
    
    integer, parameter :: nint = 1
    real(8), parameter :: pi = 3.1415926535D0
    real(8) :: DD, FF, GG, HH
    real(8) :: D(4,4), Be(4,12), Ne(2,12), S(12,12), M(12,12), UE(6), P(4)
    real(8) :: detJ, rr, ri, si
    integer :: ITYPE, IST
    integer :: JJ, lx, ly

    NPAR1  = NPAR(1)
    NUME   = NPAR(2)
    NUMMAT = NPAR(3) 
    ITYPE  = NPAR(4)

    ND = 6

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
                    4 (' .'),'( NPAR(3) ) . . =',I5,//, &
                    ' PROBLEM TYPE ',13(' .'),'( NPAR(4) ) . . =',I5,/,   &
                    '     EQ.0, AXISYMMETRIC',/,      &
                    '     EQ.1, PLANE STRAIN',/,  &
                    '     EQ.2, PLANE STRESS',/)") NUMMAT, ITYPE

        WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL    THIC',/,  &
                    ' NUMBER     MODULUS',10X,'AREA',10X,'VALUE'/,  &
                    15 X,'E',14X,'A',14X,'T')")
     
        
        DO I=1,NUMMAT
            READ (IIN,'(I5,3F10.0)') N,E(N),PR(N),THIC(N)                      ! Read material information  ÃÌº”≤¥À…±»°¢∫Ò∂»
            WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,2X,E12.5)") N,E(N),PR(N),THIC(N)
        END DO
        
        WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                        ' ELEMENT     NODE     NODE     NODE       MATERIAL',/,   &
                        ' NUMBER-N    II(1)    II(2)    II(3)     SET NUMBER')")

        N=0
        DO WHILE (N .NE. NUME)
            READ (IIN,'(5I5)') N,II,MTYPE  ! Read in element information

!           Save element information
            XY(1,N)=X(II(1))  ! Coordinates of the element's node 1
            XY(2,N)=Y(II(1))

            XY(3,N)=X(II(2))  ! Coordinates of the element's node 2
            XY(4,N)=Y(II(2))
            
            XY(5,N)=X(II(3))  ! Coordinates of the element's node 3
            XY(6,N)=Y(II(3))

            MATP(N)=MTYPE     ! Material type

            DO L=1,6
                LM(L,N)=0
            END DO

            DO L=1,2
                DO J = 1,3
                    LM(L+2*(J-1),N) = ID(L,II(J))     ! Connectivity matrix
                END DO
            END DO

!           Update column heights and bandwidth
            CALL COLHT (MHT,ND,LM(1,N))   

            WRITE (IOUT,"(I5,6X,I5,3(4X,I5),7X,I5)") N,II,MTYPE
            WRITE (10,"(3(4X,I5))") II
        END DO

        RETURN
        
    ELSE IF (IND == 2) THEN
        
        DO n = 1, nume
            mtype = MATP(N)
            
!           Obtaint stress-strain law
            FF = E(mtype) / (1. + PR(mtype))
            GG = FF * PR(mtype) / (1. - 2. * PR(mtype))
            HH = FF + GG

!           Plane strain analysis
            D(1,1)=HH    
            D(1,2)=GG
            D(1,3)=0.      
            D(2,1)=GG
            D(2,2)=HH
            D(2,3)=0.
            D(3,1)=0.
            D(3,2)=0.
            D(3,3)=FF/2. 

            IF (ITYPE .EQ. 1) THEN
                THIC(mtype) = 1.0
            END IF
            
!           Axisymmetric analysis
            D(1,4)=GG
            D(2,4)=GG
            D(3,4)=0.
            D(4,1)=GG
            D(4,2)=GG
            D(4,3)=0.
            D(4,4)=HH
            
!           For plane stress analysis condense stress-stain matrix
            IF (ITYPE == 2) THEN
                DO I = 1,3
                    DD = D(I,4) / D(4,4)
                    DO J = I,3
                        D(I,J) = D(I,J) - D(4,J) * DD
                        D(J,I) = D(I,J)
                    END DO
                END DO 
            END IF
            
!           Calculate element stiffness
            S = 0.
            M = 0.
            JJ = 0
            IF (ITYPE == 0) THEN    !Axisymmetric
                ri = (xy(1,n)+xy(3,n)+xy(5,n)) / 3
                si = (xy(2,n)+xy(4,n)+xy(6,n)) / 3

!               evaluate derivative operator B and the Jacobian determinant detJ
                call BNmat(xy(1,n), ri, si, Be, Ne, detJ, rr)
                S = S + 2 * pi * rr * &
                        matmul(transpose(Be),matmul(D,Be)) * detJ
!               The mass matrix (NOTE: NOT HAVE THE DENSITY)
                M = M + 2 * pi * rr * &
                        matmul(transpose(Ne),Ne) * detJ
            ELSE
                ri = (xy(1,n)+xy(3,n)+xy(5,n)) / 3
                si = (xy(2,n)+xy(4,n)+xy(6,n)) / 3
                    
!               evaluate derivative operator B and the Jacobian determinant detJ
                call BNmat(xy(1,n), ri, si, Be, Ne, detJ, rr)
                S = S + matmul(transpose(Be(1:3,1:12)),matmul(D(1:3,1:3),Be(1:3,1:12))) * detJ
!               The mass matrix (NOTE: NOT HAVE THE DENSITY)
                M = M + matmul(transpose(Ne),Ne) * detJ
            END IF
            
            CALL ADDBAN (DA(NP(3)),IA(NP(2)),S(1:6,1:6),LM(1,N),ND)
            CALL ADDBAN (DA(NP(13)),IA(NP(2)),M(1:6,1:6),LM(1,N),ND)
            
        END DO
        
        RETURN
        
    ELSE IF (IND == 3) THEN
        IPRINT=0
        DO N=1,NUME
            IPRINT=IPRINT + 1
            IF (IPRINT.GT.50) IPRINT=1
            IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                                'E L E M E N T  G R O U P',I4,//,   &
                                                '  ELEMENT',5X,'GAUSS POINT',5X,'StressXX',5X,'StressYY',5X,'StressXY',/,&
                                                '  NUMBER')") NG
            MTYPE=MATP(N)
            DO L = 1,3
                I = LM(2*L-1,N)
                IF (I .GT. 0) THEN
                    UE(2*L-1) = U(I)
                ELSE
                    UE(2*L-1) = 0
                END IF     
            
                J = LM(2*L,N)
                IF (J .GT. 0) THEN
                    UE(2*L) = U(J)
                ELSE
                    UE(2*L) = 0
                END IF
            END DO
        
            IF (ITYPE == 0) THEN    !Axisymmetric
                ri = (xy(1,n)+xy(3,n)+xy(5,n)) / 3
                si = (xy(2,n)+xy(4,n)+xy(6,n)) / 3
                    
!               evaluate derivative operator B and the Jacobian determinant detJ
                call BNmat(xy(1,n), ri, si, Be, Ne, detJ, rr)
                P = matmul(D,matmul(Be(:,1:6),UE))
                write (IOUT,"(I5,5X,f6.3,2X,f6.3,4X,E13.6,4X,E13.6,4X,E13.6)")N,ri,si,P(1),P(2),P(3),P(4)
            ELSE
                ri = (xy(1,n)+xy(3,n)+xy(5,n)) / 3
                si = (xy(2,n)+xy(4,n)+xy(6,n)) / 3

!               evaluate derivative operator B and the Jacobian determinant detJ
                call BNmat(xy(1,n), ri, si, Be, Ne, detJ, rr)
                P(1:3) = matmul(D(1:3,1:3),matmul(Be(1:3,1:6),UE))
                write (IOUT,"(I5,5X,f6.3,2X,f6.3,4X,E13.6,4X,E13.6,4X,E13.6)")N,ri,si,P(1),P(2),P(3)
            END IF
        END DO
    ELSE 
        STOP "*** ERROR *** Invalid IND value."
    END IF

END SUBROUTINE triangle6


subroutine BNmat6 (xy,ri,si,Be,Ne,detJ,r)
    implicit none
    
    real(8) :: xy(6), Be(4,12), Ne(2,12)
    real(8) :: ri, si
    
    real(8) :: xx(3), yy(3), xsi(3), N(6)
    real(8) :: detJ, r

    integer :: i
    
    xx(1) = xy(3) - xy(5)
    xx(2) = xy(5) - xy(1)
    xx(3) = xy(1) - xy(3)
    yy(1) = xy(4) - xy(6)
    yy(2) = xy(6) - xy(2)
    yy(3) = xy(2) - xy(4)
    
    detJ = (xy(3)*xy(6)-xy(5)*xy(4) - xy(1)*xy(6)+xy(5)*xy(2) + xy(1)*xy(4)-xy(3)*xy(2))
    
    xsi(1) = (xy(3)*xy(6) - xy(5)*xy(4) + yy(1)*ri - xx(1)*si) / detJ
    xsi(2) = (xy(5)*xy(2) - xy(1)*xy(6) + yy(2)*ri - xx(2)*si) / detJ
    xsi(3) = (xy(1)*xy(4) - xy(3)*xy(2) + yy(3)*ri - xx(3)*si) / detJ
    
    N(1) = xsi(1)*(2*xsi(1)-1)
    N(2) = xsi(2)*(2*xsi(2)-1)
    N(3) = xsi(3)*(2*xsi(3)-1)
    N(4) = 4*xsi(1)*xsi(2)
    N(5) = 4*xsi(2)*xsi(3)
    N(6) = 4*xsi(3)*xsi(1)
    
    Ne = 0

    do i = 1,6
        Ne(1,2*i-1) = N(i)
        Ne(2,2*i-1) = 0
        Ne(1,2*i  ) = 0
        Ne(2,2*i  ) = N(i)
    end do
    
    r = xsi(1)*xy(1) + xsi(2)*xy(3) + xsi(3)*xy(5)
    
    Be = 0
    
    do i = 1,3
        Be(1,2*i-1) =  yy(i)*(4*xsi(i)-1)/detJ 
        Be(3,2*i-1) = -xx(i)*(4*xsi(i)-1)/detJ
        Be(4,2*i-1) =  N(i)/r
    
        Be(2,2*i) = -xx(i)*(4*xsi(i)-1)/detJ
        Be(3,2*i) =  yy(i)*(4*xsi(i)-1)/detJ
    end do
    
    Be(1,7) =  4*(xsi(2)*yy(1) + xsi(1)*yy(2)) / detJ
    Be(3,7) = -4*(xsi(2)*xx(1) + xsi(1)*xx(2)) / detJ
    Be(4,7) =  N(4) / r
    Be(2,8) = -4*(xsi(2)*xx(1) + xsi(1)*xx(2)) / detJ
    Be(3,8) =  4*(xsi(2)*yy(1) + xsi(1)*yy(2)) / detJ
    
    Be(1,9) =  4*(xsi(2)*yy(3) + xsi(3)*yy(2)) / detJ
    Be(3,9) = -4*(xsi(2)*xx(3) + xsi(3)*xx(2)) / detJ
    Be(4,9) =  N(5) / r
    Be(2,10) = -4*(xsi(2)*xx(3) + xsi(3)*xx(2)) / detJ
    Be(3,10) =  4*(xsi(2)*yy(3) + xsi(3)*yy(2)) / detJ

    Be(1,11) =  4*(xsi(1)*yy(3) + xsi(3)*yy(1)) / detJ
    Be(3,11) = -4*(xsi(1)*xx(3) + xsi(3)*xx(1)) / detJ
    Be(4,11) =  N(6) / r
    Be(2,12) = -4*(xsi(1)*xx(3) + xsi(3)*xx(1)) / detJ
    Be(3,12) =  4*(xsi(1)*yy(3) + xsi(3)*yy(1)) / detJ
    
end subroutine BNmat6