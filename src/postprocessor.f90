SUBROUTINE POSTPROCESSOR
  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TECPLOT                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  IMPLICIT NONE
  CHARACTER*80 FileInp
  INTEGER::NPAR1,ID(NDF,NUMNP),N,I,KK,II,IL,P1,P2,P3,P4,P5,P6,P7,P8,NUME
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),DX(NUMNP),DY(NUMNP),DZ(NUMNP),SXX(NUMNP),SYY(NUMNP),SXY(NUMNP),DISP(NEQ)
 
   OPEN(ITEC  , FILE = "TEC.DAT", STATUS = "REPLACE")
 
   
    NPAR1=NPAR(1)
   IF (NPAR1==2) THEN
  WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"","" D-X"" ,"" D-Y"","" D-Z"",""SXX"",""SYY"",""SXY""'//,&
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
      
       ! WRITE (ITEC,"(',J=3',I5)",advance="no") 
    
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
  
   
   Do N=1,NUMNP
     read(10,'(3E18.6)')  DX(N),DY(N),DZ(N)
     
   enddo
   
   
   IF (NPAR1==2) THEN
    DO N=1,NUMNP
    READ(10,'(3F15.6)') SXX(N),SYY(N),SXY(N)
    WRITE(ITEC,'(9E14.6)')X(N),Y(N),Z(N),DX(N),DY(N),DZ(N),SXX(N),SYY(N),SXY(N)
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
      
  ELSE
  END IF
   END DO
   
    END SUBROUTINE
      