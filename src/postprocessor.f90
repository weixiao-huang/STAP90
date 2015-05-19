SUBROUTINE POSTPROCESSOR
  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TECPLOT                                            .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  IMPLICIT NONE
  CHARACTER*80 FileInp
  INTEGER::NPAR1,ID(NDF,NUMNP),N,I,KK,II,IL,P1,P2,P3,P4,P5,P6,P7,P8,NUME
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),D(3),DISP(NEQ)
 
   OPEN(ITEC  , FILE = "TEC.DAT", STATUS = "REPLACE")
 
  WRITE (ITEC,"('TITLE="" ',A80,'""'/,&
  'VARIABLES=""X"",""Y"","" Z"","" D-X"" ,"" D-Y"","" D-Z""'//,&
  'ZONE  F=FEPOINT ')",advance="no") HED
  
  WRITE (ITEC,"(',N=',I5)",advance="no") NUMNP
      
  WRITE (ITEC,"(',E=',I5)",advance="no") NPAR(2) 
  
  NPAR1=NPAR(1)
 
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
     read(10,'(3E18.6)')  D(1),D(2),D(3)
     write(ITEC,'(6E14.6)')X(N),Y(N),Z(N),D(1),D(2),D(3)
   enddo
   
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
      