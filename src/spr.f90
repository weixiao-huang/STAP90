subroutine ADDSPR (NOE,NPELM,TRUESTR,CIGMA)

  USE GLOBALS
  USE MEMALLOCATE
  
  IMPLICIT NONE
  REAL(8) :: CIGMA(NDF,NIE+1),TRUESTR(NDF,NUMNP+NPAR(2))    ! NPAR(2)=NUME
  INTEGER :: NPELM(NIE)
  INTEGER :: N ,I, J, JJ,NOE            ! NOEΪ element NO.
  REAL(8) :: ANS(NDF,NIE+1)
  REAL(8) :: M(5,5)=reshape((/1.8660254038,-0.5,0.1339745962 ,-0.5,0.0,-0.5,1.8660254038,-0.5,0.1339745962 ,0.0,0.1339745962 ,-0.5,1.8660254038, -0.5,0.0,-0.5,0.1339745962 ,-0.5,1.8660254038,0.0,0.0,0.0,0.0,0.0,1.0/),(/5,5/))
  
  
  ANS=matmul(CIGMA,M)                   ! 单元应力磨平
  
  DO I=1,NDF
      DO J=1,NIE
           JJ=NPELM(J)
        TRUESTR(I,JJ)=TRUESTR(I,JJ)+ANS(I,J)  
      END DO

   TRUESTR(I,NUMNP+NOE)=ANS(I,NIE+1)
      
  END DO

  RETURN
  
    END SUBROUTINE ADDSPR
    
    
subroutine SPR_4Q( NPELM,TRUESTR,XYZ)

  USE GLOBALS
  USE MEMALLOCATE
  
  IMPLICIT NONE
  REAL(8) :: TRUESTR(NDF,NUMNP+NPAR(2)) ,XYZ(8,NPAR(2)), XYZGP(2,NPAR(2))       ! NPAR(2)=NUME
  INTEGER :: NPELM(NIE,NPAR(2)),NPIWE(NIE+1,NUMNP)              ! NPIWE : Nodel Point in Which element?
  INTEGER :: N ,I, J,JJ,NOE            ! NOEΪ element NO.
  REAL(8) :: ANS(NDF,NUMNP)
  
  DO I=1,NIE+1
      DO J=1,NUMNP
          NPIWE(I,J)=0.0
      enddo
  enddo
  
  Do I=1,NPAR(2)
      XYZGP(1,I)=(XYZ(1,I)+XYZ(3,I)+XYZ(5,I)+XYZ(7,I))/4.0;
      XYZGP(2,I)=(XYZ(2,I)+XYZ(4,I)+XYZ(6,I)+XYZ(8,I))/4.0;
  ENDDO
      
  
  Do I=1,NPAR(2)
      DO J=1,NIE
          JJ=NPELM(I,J)
          NPIWE(1,JJ)=NPIWE(1,JJ)+1
          NPIWE(NPIWE(1,JJ)+1,JJ)=I
      ENDDO
  ENDDO
  
   DO I=1,NDF
       DO J=1,NUMNP
           ANS(I,J)=TRUESTR(I,J)/NPIWE(1,J)
       ENDDO
   ENDDO
          
 DO I=1,NUMNP
     
    if(NPIWE(1,I)==4) then
        
        
        
        
    endif    
enddo
   
    WRITE(10,'(<NUMNP>F15.6)') ANS
  
  end subroutine SPR_4Q