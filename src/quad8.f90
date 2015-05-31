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

SUBROUTINE Quad8
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
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
      MM = 3*NUMMAT*ITWO +17*NUME + 16*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: PR(NUMMAT)
! N103: LM(16,NUME)
! N104: XYZ(16,NUME)
! N105: MTAP(NUME)
! N106: THICK(NUMMAT)

  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+16*NUME
  N105=N104+16*NUME*ITWO
  N106=N105+NUME
  N107=N106+NUMMAT*ITWO
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL QUAD8_SUB (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N106),A(N103),A(N104),A(N105))


  RETURN

END SUBROUTINE Quad8


SUBROUTINE QUAD8_SUB (ID,X,Y,Z,U,MHT,E,PR,THICK,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   8 nodes quadr element subroutine                                .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(3,NUMNP),LM(16,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),PR(NPAR(3)),  &
            THICK(NPAR(3)), XYZ(16,NPAR(2)),U(NEQ),UE(16)                           
  REAL(8) :: S(16,16),M(8,8),D(4,4),E0,PR0,D0,r
  INTEGER :: NPAR1, NUME, NUMMAT, ND,P(8), L, N, I,J,K
  INTEGER :: MTYPE, IPRINT, ITYPE                                        
  REAL(8) ::  XM, XX, YY, STR(4), PF(4)
  REAL(8) :: GP(2), WGT(2),B(4,16),detJ,NL(2,8),NQ(1,8)
  REAL(8),parameter:: pi=3.141592654
  !GP =[-0.9061798459, -0.5384693101 ,0.0 ,0.5384693101,0.9061798459]                          !五点GAUSS积分
  !WGT=[ 0.2369268851,  0.4786286705 ,0.5688888889, 0.4786286705,0.2369268851]
  GP=[-0.5773502692 , 0.5773502692]
  WGT=[ 1          ,          1]
  
  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 
  ITYPE  = NPAR(4)
  ND=16

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
                    '     EQ.8, 8Q ELEMENTS',//,      &
                  ' SOLVE MODE ',13(' .'),'( NPAR(4) ) . . =',I5,/,   &
                    '     EQ.0, AXIS SYMMETRY',/,      &            
                    '     EQ.1, PLAIN STRAIN',/,  &
                    '     EQ.2, PLAIN STRSS',//,   &    
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,ITYPE,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND PROPETIES ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S       POISSION      THICKNESS',/,  &
                   ' NUMBER     MODULUS         RATIO         VALUE  ',/,  &
                   '               E             PR             T    ')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,3F12.5)') N,E(N),PR(N),THICK(N)                      ! Read material information  
        WRITE (IOUT,"(I5,4X,E12.5,2X,E12.5,2X,E12.5)") N,E(N),PR(N),THICK(N)
     END DO
           
     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE    NODE    NODE    NODE    NODE    NODE    NODE    NODE    MATERIAL',/,   &
                      ' NUMBER-N     P1      P2      P3      P4      P5      P6      P7      P8     SET NUMBER')")
     
     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(10I5)') N,P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8),MTYPE  ! Read in element information

!       Save element information
        XYZ(1,N)=X(P(1))     ! Coordinates of the element's first node
        XYZ(2,N)=Y(P(1))
        XYZ(3,N)=X(P(2))     ! Coordinates of the element's second node
        XYZ(4,N)=Y(P(2))      
        XYZ(5,N)=X(P(3))     ! Coordinates of the element's third node
        XYZ(6,N)=Y(P(3))
        XYZ(7,N)=X(P(4))     ! Coordinates of the element's fourth node
        XYZ(8,N)=Y(P(4))
        XYZ(9,N)=X(P(5))     ! Coordinates of the element's fifth node
        XYZ(10,N)=Y(P(5))
        XYZ(11,N)=X(P(6))     ! Coordinates of the element's sixth node
        XYZ(12,N)=Y(P(6))      
        XYZ(13,N)=X(P(7))     ! Coordinates of the element's seventh node
        XYZ(14,N)=Y(P(7))
        XYZ(15,N)=X(P(8))     ! Coordinates of the element's eighth node
        XYZ(16,N)=Y(P(8))

        MATP(N)=MTYPE  ! Material type
        
        DO L=1,16
           LM(L,N)=0
        END DO
              
        DO L=1,2
           LM(L,N)=ID(L,P(1))     ! Connectivity matrix
           LM(L+2,N)=ID(L,P(2)) 
           LM(L+4,N)=ID(L,P(3)) 
           LM(L+6,N)=ID(L,P(4))
           LM(L+8,N)=ID(L,P(5))
           LM(L+10,N)=ID(L,P(6)) 
           LM(L+12,N)=ID(L,P(7)) 
           LM(L+14,N)=ID(L,P(8))

        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,7X,I5)") N,P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8),MTYPE
        WRITE (10,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)") P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8)
     END DO
     
     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN


     DO N=1,NUME
        MTYPE=MATP(N)
        D = 0
        IF(ITYPE == 1 .OR. ITYPE == 0) THEN   
            E0 = E(MTYPE)/(1-PR(MTYPE)*PR(MTYPE))     !!!!!!!!!!!!!!平面应变
            PR0 = PR(MTYPE)/(1-PR(MTYPE))
            THICK(MTYPE) = 1.
            IF(ITYPE == 0) THEN
                D(1,4) = PR0
                D(2,4) = PR0
                D(4,1) = D(1,4)
                D(4,2) = D(2,4)
                D(4,4) = 1
            END IF
        ELSE IF(ITYPE == 2) THEN
            E0 = E(MTYPE)                               !!!!!!!!!!!!平面应力
            PR0 = PR(MTYPE)
        ENDIF
        D0 = E0/(1-PR0*PR0)
        D(1,1) = 1
        D(1,2) = PR0
        D(2,1) = D(1,2)
        D(2,2) = 1
        D(3,3) = 0.5*(1-PR0)
        D = D0*D
            
        S = 0
   
        if(ITYPE.eq.1.or.ITYPE.eq.2) then
            do i=1,2                                                      
                do j=1,2
                    CAll Bmat_8Q(gp(i),gp(j),XYZ(1,N),B,detJ)
                    CALL Nmat_8Q(gp(i),gp(j),NQ)
                    S=S+WGT(i)*WGT(j)*matmul(matmul(transpose(B),D),B)*detJ*THICK(MTYPE)
                    
!                   The mass matrix (NOTE: NOT HAVE THE DENSITY)
                    M=M+WGT(i)*WGT(j)*matmul(transpose(NQ),NQ)*detJ*THICK(MTYPE)*1e3
                end do
            end do
        else if(ITYPE.eq.0)then                                                  !轴对称问题
            do i=1,2                                                      
                do j=1,2
                    CAll Bmat_8Q(gp(i),gp(j),XYZ(1,N),B,detJ)
                    CALL Nmat_8Q(gp(i),gp(j),NQ)                    
                    CALL Nmat(gp(i),gp(j),NL)
                    r=2*pi*(NL(1,1)*XYZ(1,N)+NL(1,3)*XYZ(3,N)+NL(1,5)*XYZ(5,N)+NL(1,7)*XYZ(7,N))
                    S=S+r*WGT(i)*WGT(j)*matmul(matmul(transpose(B),D),B)*detJ
                    
!                   The mass matrix (NOTE: NOT HAVE THE DENSITY)
                    M=M+r*WGT(i)*WGT(j)*matmul(transpose(NQ),NQ)*detJ*1e3
                end do
            end do
        end if       

        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)
        CALL ADDBAN (DA(NP(13)),IA(NP(2)),M,LM(1,N),ND)
     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     CALL MEMALLOC(9,"NPFORCE",(NUME+NUMNP)*NDF,ITWO)                                   !分配内存9 以储存节点及用于重构的单元超收敛点的位移

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                            'E L E M E N T  G R O U P',I4,//,   &
                                            '  ELEMENT',5X,'GAUSS POINT',5X,'StressXX',5X,'StressYY',5X,'StressXY',/,&
                                            '  NUMBER')") NG
        MTYPE=MATP(N)
        DO L=1,8    
            I=LM(2*L-1,N)
            if (I.GT.0)then
                UE(2*L-1)=U(I)
            else
                UE(2*L-1)=0
            endif       
            J=LM(2*L,N)
            IF (J.GT.0)then
                UE(2*L)=U(J)
            else
                UE(2*L)=0
            endif
        END DO
        
        do i=1,2                                                      
            do j=1,2  
                CAll Bmat_8Q(gp(i),gp(j),XYZ(1,N),B,detJ)
                STR = matmul(B,UE)
                PF  = matmul(D,STR)
                WRITE (IOUT,"(I5,5X,f6.3,2X,f6.3,4X,E13.6,4X,E13.6,4X,E13.6)")N,gp(i),gp(j),PF(1),PF(2),PF(3)
            end do
        end do
 
  
     END DO

  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

END SUBROUTINE QUAD8_SUB


subroutine Nmat_8Q(xi,eta,N)
    implicit none
    real(8) :: xi,eta,N(8)
    N(1)=-0.25*(1-xi)*(1-eta)*(xi+eta+1)
    N(2)=0.25*(1+xi)*(1-eta)*(xi-eta-1)
    N(3)=0.25*(1+xi)*(1+eta)*(xi+eta-1)
    N(4)=-0.25*(1-xi)*(1+eta)*(xi-eta+1)
    N(5)=0.5*(1+xi)*(1-xi)*(1-eta)
    N(6)=0.5*(1+xi)*(1+eta)*(1-eta)
    N(7)=0.5*(1+xi)*(1-xi)*(1+eta)
    N(8)=0.5*(1-xi)*(1+eta)*(1-eta)
end subroutine Nmat_8Q
    
subroutine Bmat_8Q(xi,eta,XY,B,detJ)
    implicit none
    real(8) :: xi,eta,XY(16),loca(8,2),B(4,16),detJ,GNxe(2,8),GNxy(2,8),J(2,2),JINV(2,2),DUM,N(8),NL(2,8)
    integer:: K,I

    GNxe(1,1)=0.25*(1-eta)*(2*xi+eta)
    GNxe(1,2)=0.25*(1-eta)*(2*xi-eta)
    GNxe(1,3)=0.25*(1+eta)*(2*xi+eta)
    GNxe(1,4)=0.25*(1+eta)*(2*xi-eta)
    GNxe(1,5)=-xi*(1-eta)
    GNxe(1,6)=0.5*(1-eta**2)
    GNxe(1,7)=-xi*(1+eta)
    GNxe(1,8)=-0.5*(1-eta**2)
    GNxe(2,1)=0.25*(1-xi)*(xi+2*eta)
    GNxe(2,2)=0.25*(1+xi)*(-xi+2*eta)
    GNxe(2,3)=0.25*(1+xi)*(xi+2*eta)
    GNxe(2,4)=0.25*(1-xi)*(-xi+2*eta)
    GNxe(2,5)=-0.5*(1-xi**2)
    GNxe(2,6)=-eta*(1+xi)
    GNxe(2,7)=0.5*(1-xi**2)
    GNxe(2,8)=-eta*(1-xi)

    DO K=0,7
        loca(K+1,1)=XY(2*K+1)
        loca(K+1,2)=XY(2*K+2)
    END DO
        
    J=matmul(GNxe,loca)
    detJ=J(1,1)*J(2,2)-J(2,1)*J(1,2)
    DUM=1./detJ
    JINV(1,1)=J(2,2)*DUM
    JINV(1,2)=-J(1,2)*DUM
    JINV(2,1)=-J(2,1)*DUM
    JINV(2,2)=J(1,1)*DUM

    GNxy=matmul(JINV,GNxe)
    CALL Nmat_8Q(xi,eta,N)
    CALL Nmat(xi,eta,NL)
    
    B=0
    DO K=0,7
        B(1,2*K+1)=GNxy(1,K+1)
        B(2,2*K+2)=GNxy(2,K+1)
        B(3,2*K+1)=GNxy(2,K+1)
        B(3,2*K+2)=GNxy(1,K+1)
        B(4,2*K+2)=N(K+1)/(NL(1,1)*XY(1)+NL(1,3)*XY(3)+NL(1,5)*XY(5)+NL(1,7)*XY(7))
    END DO

end subroutine Bmat_8Q