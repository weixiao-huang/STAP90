SUBROUTINE POSTPROCESSOR
  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TECPLOT                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  IMPLICIT NONE
  CHARACTER*80 FileInp
  INTEGER::NPAR1,ID(NDF,NUMNP),N,I,KK,II,IL,P1,P2,P3,P4,P5,P6,P7,P8,NUME
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),DX(NUMNP),DY(NUMNP),DZ(NUMNP),SXX(NUMNP),SYY(NUMNP),SXY(NUMNP),DISP(NEQ),CX(NUMNP),CY(NUMNP)
 
   OPEN(ITEC  , FILE = "TEC.DAT", STATUS = "REPLACE")
 
   
    NPAR1=NPAR(1)
   IF (NPAR1==2) THEN
  WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"","" D-X"" ,"" D-Y"","" D-Z"",""SXX"",""SYY"",""SXY""'//,&
  'ZONE  F=FEPOINT ')",advance="no") HED
  ELSE IF(NPAR1==6) THEN
  WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"","" W"" ,"" CX"","" CY""'//,&
  'ZONE  F=FEPOINT ')",advance="no") HED
   ELSE IF(NPAR1==7) THEN
  WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"",""D-X"" ,"" D-Y"","" D-Z"",""CX"",""CY""'//,&
  'ZONE  F=FEPOINT ')",advance="no") HED
  
    ELSE
   WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"","" D-X"" ,"" D-Y"","" D-Z""'//,&
  'ZONE  F=FEPOINT ')",advance="no") HED
   END IF
        
        
  WRITE (ITEC,"(',N=',I5)",advance="no") NUMNP
      
  WRITE (ITEC,"(',E=',I5)",advance="no") NPAR(2) 
  
 
 
  IF (NPAR1 == 2) THEN    ! Quadrilateral Elements
      WRITE (ITEC,"(',ET=QUADRILATERAL')",advance="no")
  ELSE IF (NPAR1 == 3) THEN    ! Triangle Elements
     WRITE (ITEC,"(',ET=TRIANGLE')",advance="no")
  ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
       WRITE (ITEC,"(',ET=BRICK')",advance="no")
  ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       REWIND(ITEC)
        WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
              'VARIABLES=""X"",""Y"","" Z"","" D-X"" ,"" D-Y"","" D-Z""'//,&
              'ZONE  F=POINT ')",advance="no") HED
  
        WRITE (ITEC,"(',I=',I5)",advance="no") NUMNP
      
  ELSE IF (NPAR1 == 6) THEN    ! Quadrilateral Elements
      WRITE (ITEC,"(',ET=QUADRILATERAL')",advance="no")  
  ELSE IF (NPAR1 == 7) THEN    ! Quadrilateral Elements
      WRITE (ITEC,"(',ET=QUADRILATERAL')",advance="no")
  ELSE
  END IF

  WRITE (ITEC,"(' C=CYAN',/)")

  
REWIND(10)

   DO N=1,NUMNP
    READ(10,"(3F13.3)") X(N),Y(N),Z(N)
   END DO
   
   DO N=1,NPAR(2)
       read (10,'()')
   enddo
 
   
   rewind(ITEC)
      DO N=1,5
       read (ITEC,'()')
   end do
  
  IF(NPAR1==7)THEN 
   Do N=1,NUMNP
     read(10,'(5E18.6)') DX(N),DY(N),DZ(N),CX(N),CY(N)
   enddo
  ELSE IF(NPAR1==6)THEN 
   Do N=1,NUMNP
     read(10,'(5E18.6)') DX(N),DY(N),DZ(N),CX(N),CY(N)
   enddo
   ELSE
   Do N=1,NUMNP
     read(10,'(3E18.6)') DX(N),DY(N),DZ(N)
     enddo
   END IF
   
   
   IF (NPAR1==2) THEN
    DO N=1,NUMNP
    READ(10,'(3F15.6)') SXX(N),SYY(N),SXY(N)
    WRITE(ITEC,'(9E14.6)')X(N),Y(N),Z(N),DX(N),DY(N),DZ(N),SXX(N),SYY(N),SXY(N)
    END DO
   ELSE IF(NPAR1==7) THEN
    DO N=1,NUMNP
    WRITE(ITEC,'(8E14.6)')X(N),Y(N),Z(N),DX(N),DY(N),DZ(N),CX(N),CY(N)
    END DO    
     ELSE IF(NPAR1==6) THEN
    DO N=1,NUMNP
    WRITE(ITEC,'(6E14.6)')X(N),Y(N),Z(N),DZ(N),CX(N),CY(N)
    END DO 
   ELSE    
       DO N=1,NUMNP
       WRITE(ITEC,'(6E14.6)')X(N),Y(N),Z(N),DX(N),DY(N),DZ(N)
       END DO
   ENDIF
   
   
   
   
    
   REWIND(10)
  
   
   DO N=1,NUMNP
    read (10,'()')
   enddo
   
   
 DO N=1,NPAR(2)
   IF (NPAR1 == 2) THEN 
     
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)")P1,P2,P3,P4
       WRITE(ITEC,"(I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4
      
   ELSE IF (NPAR1 == 3) THEN! Triangle Elements
    
       READ(10,"(3(4X,I5))")P1,P2,P3
       WRITE(ITEC,"(I5,3X,I5,3X,I5)")P1,P2,P3
    ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
        
       READ(10,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
       WRITE(ITEC,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
    ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       
       READ(10,"(I5,4X,I5)") P1,P2
       WRITE(ITEC,"(I5,4X,I5,4X,I5)") P1,P2
    ELSE IF (NPAR1 == 6) THEN    ! PLATE Elements
       
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4 
       WRITE(ITEC,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4 

     ELSE IF (NPAR1 == 7) THEN    ! shell Elements
       
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4 
       WRITE(ITEC,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4
       
    ELSE IF (NPAR1 == 9) THEN    ! 8Q Elements
       
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4 
       WRITE(ITEC,"(I5,4X,I5,4X,I5,4X,I5)") P1,P2,P3,P4
      
  ELSE
  END IF
   END DO
   
    END SUBROUTINE
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    

SUBROUTINE POSTEIGEN
  USE GLOBALS
  IMPLICIT NONE
  CHARACTER*80 FileInp
  INTEGER::NPAR1,ID(NDF,NUMNP),N,I,KK,II,IL,P1,P2,P3,P4,P5,P6,P7,P8,NUME,P,VV
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),D(3),DISP(NEQ),R(3),DY(3,NUMNP),RY(3,NUMNP),M(3,NUMNP),MX(3)
 
  OPEN(16,FILE='EIGEN.DAT',STATUS='REPLACE')
  
  REWIND(10)
 


  
  NPAR1=NPAR(1)
  
  WRITE (16,"('TITLE="" ',A80,'""'/,&
              'VARIABLES=""X"",""Y"","" Z""'//,&
              'ZONE  F=FEPOINT ')",advance="no") HED
  
  WRITE (16,"(',N=',I5)",advance="no") NUMNP
  WRITE (16,"(',E=',I5)",advance="no")NPAR(2)
  IF (NPAR1 == 2) THEN    ! Quadrilateral Elements
   WRITE (16,"(',ET=QUADRILATERAL')")
    ELSE IF (NPAR1 == 3) THEN    ! Triangle Elements
     WRITE (16,"(',ET=TRIANGLE')")
    ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
       WRITE (16,"(',ET=BRICK')")
    ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       REWIND(16)
        WRITE (16,"('TITLE="" ',A80,'""'/,&
              'VARIABLES=""X"",""Y"","" Z""'//,&
              'ZONE  F=POINT ')",advance="no") HED 
        WRITE (16,"(',I=',I5)") NUMNP
    END IF

     DO N=1,NUMNP
       READ(10,'(3E13.3)')D(1),D(2),D(3)
       WRITE(16,'(3E13.3)')D(1),D(2),D(3)
        DO I=1,3
         DY(I,N)=D(I)
        END DO
     END DO   
     
     
    DO N=1,NPAR(2)
   IF (NPAR1 == 2) THEN 
     
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)")P1,P2,P3,P4
       WRITE(16,"(I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4
      
   ELSE IF (NPAR1 == 3) THEN! Triangle Elements
    
       READ(10,"(3(4X,I5))")P1,P2,P3
       WRITE(16,"(I5,3X,I5,3X,I5)")P1,P2,P3
    ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
        
       READ(10,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
       WRITE(16,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
    ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       
       READ(10,"(I5,4X,I5)") P1,P2
       WRITE(16,"(I5,4X,I5,4X,I5)") P1,P2
      
  ELSE
  END IF
 END DO  
       
 DO P=1,NPAR(6)   
     REWIND(10)
     DO I=1,3
         DO N=1,NUMNP
        RY(I,N)=0
         END DO
     END DO
     
  IF(NPAR1==2)THEN
      DO  N=1,(P+2)*NUMNP+NPAR(2) 
        READ(10,"()")
      END DO
  ELSE
      DO  N=1,(P+1)*NUMNP+NPAR(2) 
      READ(10,"()")
      END DO
   END IF
  
  DO N=1,NUMNP
      READ(10,"(3F18.6)")R(1),R(2),R(3)
       DO I=1,3
         RY(I,N)=R(I)
       END DO
  END DO   
   
  M=DY+RY
  
  DO N=1,(NPAR(6)-1)*NUMNP
      READ(10,'()')
  END DO

  WRITE(10,"(3F18.6)")M
END DO
 
 DO P=1,NPAR(6)
      WRITE (16,"('ZONE C=CYAN,F=FEPOINT')",advance="no") 
      WRITE (16,"(',N=',I5)",advance="no") NUMNP
      
      WRITE (16,"(',E=',I5)",advance="no") NPAR(2) 

    IF (NPAR1 == 2) THEN    ! Quadrilateral Elements
      WRITE (16,"(',ET=QUADRILATERAL')")
    ELSE IF (NPAR1 == 3) THEN    ! Triangle Elements
     WRITE (16,"(',ET=TRIANGLE')")
    ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
       WRITE (16,"(',ET=BRICK')")
    ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       REWIND(16)
        WRITE (16,"('TITLE="" ',A80,'""'/,&
              'VARIABLES=""X"",""Y"","" Z""'//,&
              'ZONE  F=POINT ')",advance="no") HED
  
        WRITE (16,"(',I=',I5)") NUMNP
    endif

   
    REWIND(10)
    IF(NPAR(1)==2)THEN
    VV=NPAR(2)+(NPAR(6)+2+P)*NUMNP
    ELSE
     VV=NPAR(2)+(NPAR(6)+1+P)*NUMNP
     END IF
        
    DO N=1,VV
        read(10,'()')
    END DO
    
   DO N=1,NUMNP
    read (10,'(3F18.6)') MX(1),MX(2),MX(3)
    WRITE(16,'(3F16.8)') MX(1),MX(2),MX(3)
   enddo
   
 REWIND(10)
 DO N=1,NUMNP
     READ(10,'()')
 END DO
 
 DO N=1,NPAR(2)
   IF (NPAR1 == 2) THEN 
     
       READ(10,"(I5,4X,I5,4X,I5,4X,I5)")P1,P2,P3,P4
       WRITE(16,"(I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4
      
   ELSE IF (NPAR1 == 3) THEN! Triangle Elements
    
       READ(10,"(3(4X,I5))")P1,P2,P3
       WRITE(16,"(I5,3X,I5,3X,I5)")P1,P2,P3
    ELSE IF (NPAR1 == 4) THEN    ! 8H Elements
        
       READ(10,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
       WRITE(16,"(I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5,3X,I5)")P1,P2,P3,P4,P5,P6,P7,P8
    ELSE IF (NPAR1 == 5) THEN    ! BEAM Elements
       
       READ(10,"(I5,4X,I5)") P1,P2
       WRITE(16,"(I5,4X,I5,4X,I5)") P1,P2
      
  ELSE
  END IF
 END DO
    
 END DO


  
  
  
   END SUBROUTINE
