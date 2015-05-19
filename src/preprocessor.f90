!---------------------------------------------------------------
!           将ANSYS的输出文件转化成STAP90.IN文件
!  已完成3T/4Q单元,ANSYS输入文件见ANS.IN,输出文件为STAP90.IN
!---------------------------------------------------------------	
SUBROUTINE PREOPENFILES()
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   采用前处理时的打开文件的子函数                                  .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS

  IMPLICIT NONE
  LOGICAL :: EX
  CHARACTER*80 FileInp
  
  CALL PREPROCESSOR

  INQUIRE(FILE = 'STAP90.IN', EXIST = EX)
  IF (.NOT. EX) THEN
     PRINT *, "*** STOP *** FILE STAP90.IN DOES NOT EXIST !"
     STOP
  END IF

  OPEN(IIN   , FILE = "STAP90.IN",  STATUS = "OLD")
  OPEN(IOUT  , FILE = "STAP90.OUT", STATUS = "REPLACE")
  OPEN(10  , FILE = "DATA.OUT", STATUS = "REPLACE")
  OPEN(IELMNT, FILE = "ELMNT.TMP",  FORM = "UNFORMATTED")
  OPEN(ILOAD , FILE = "LOAD.TMP",   FORM = "UNFORMATTED")
  
    END SUBROUTINE PREOPENFILES
    
    
SUBROUTINE PREPROCESSOR
!---------------------------------------------------------------
!           将ANSYS的输出文件转化成STAP90.IN文件
!  已完成3T/4Q单元,ANSYS输入文件见ANS.IN,输出文件为STAP90.IN
!---------------------------------------------------------------	
IMPLICIT NONE

INTEGER :: NUMP,NUME,NUMMAT,NUMEG
REAL(8),ALLOCATABLE :: POS(:,:)       !节点坐标矩阵
INTEGER,ALLOCATABLE :: IEN(:,:)       !单元IEN矩阵
INTEGER :: ELETY                      !单元类型
INTEGER,ALLOCATABLE :: ID(:,:)        !各个自由度是否被约束
REAL(8),ALLOCATABLE :: F(:,:)         !各个节点各方向是否有集中载荷
INTEGER :: NPAR1
TYPE MATERIAL
    REAL(8)::E
    REAL(8)::P
    REAL(8)::DENS
    REAL(8)::AREA
END TYPE
TYPE(MATERIAL)::MAT
INTEGER::I,J
INTEGER::NUMC
REAL(8)::FORCE
CHARACTER(6)::P,P1,P2
INTEGER::IANS=1
INTEGER::IIN=2
OPEN(IIN,FILE="STAP90.IN",STATUS="UNKNOWN")
OPEN(IANS,FILE="ANS.IN",STATUS="OLD")
!读入单元，节点，材料，单元组个数
DO WHILE (.TRUE.)
	READ (IANS,'(A5)') P
	IF (P=='*ELSE') THEN
		READ (IANS,'(12X,I9)') NUMP
		READ (IANS,'(12X,I9)') NUME
		EXIT
	ELSE
		CYCLE
	ENDIF
ENDDO
!读入单元类型
DO WHILE (.TRUE.)
	READ (IANS,'(A3)') P
	IF (P=='DOF') THEN
	    READ (IANS,'(12X,I3)') ELETY
	    EXIT
	ELSE
	    CYCLE
	ENDIF		
ENDDO
!读入各个节点的坐标
ALLOCATE(POS(3,NUMP))
DO WHILE (.TRUE.)
	READ (IANS,'(A6)') P
	IF (P=='(3i9,6') THEN
		DO I=1,NUMP
			READ (IANS,'(27X,3F20.13)') POS(1,I),POS(2,I),POS(3,I)
		ENDDO
		EXIT
	ELSE 
	    CYCLE
	ENDIF
ENDDO
!读入单元的节点的全局编号
IF (ELETY==180) THEN
    ALLOCATE(IEN(2,NUME))
ELSEIF (ELETY==182) THEN
    ALLOCATE(IEN(4,NUME))
ENDIF
DO WHILE (.TRUE.)
	READ (IANS,'(A6)') P
	IF (P=='(19i9)') THEN
	    DO I=1,NUME
	        IF (ELETY==182) THEN
 			    READ (IANS,'(99X,4I9)') IEN(1,I),IEN(2,I),IEN(3,I),IEN(4,I)
 			ELSEIF (ELETY==180) THEN
 			    READ (IANS,'(99X,2I9)') IEN(1,I),IEN(2,I)
 			ENDIF
		ENDDO
		EXIT
	ELSE
		CYCLE
	ENDIF
ENDDO
!读入边界条件
ALLOCATE(ID(3,NUMP),F(3,NUMP))
DO I=1,NUMP
    ID(1,I)=0;ID(2,I)=0;ID(3,I)=1
    F=0
ENDDO
NUMC=0    !计数集中载荷个数
DO WHILE (.TRUE.)
	READ (IANS,'(A5)') P
	IF (P=='ERESX') THEN
	    !读取第一个要先换一行
	    READ (IANS,'(/,A1)',ADVANCE='NO') P2
		IF (P2=='D') THEN
			READ (IANS,'(X,I7,X,A2)') I,P1
			IF (P1=='UX') THEN
				ID(1,I)=1
			ELSEIF (P1=='UY') THEN
			    ID(2,I)=1
			ELSEIF (P1=='UZ') THEN
				ID(3,I)=1
			ENDIF
		ELSEIF (P2=='F') THEN
			READ (IANS,'(X,I7,X,A2,3X,F12.4)') I,P1,FORCE
			IF (FORCE/=0) NUMC=NUMC+1
			IF (P1=='FX') THEN
				F(1,I)=FORCE
			ELSEIF (P1=='FY') THEN
				F(2,I)=FORCE
			ELSEIF (P1=='FZ') THEN
				F(3,I)=FORCE
			ENDIF
		ENDIF
		DO WHILE (.TRUE.)
			READ (IANS,'(A1)',ADVANCE='NO') P2
			IF (P2=='D') THEN
				READ (IANS,'(X,I7,X,A2)') I,P1
				IF (P1=='UX') THEN
					ID(1,I)=1
				ELSEIF (P1=='UY') THEN
					ID(2,I)=1
				ELSEIF (P1=='UZ') THEN
					ID(3,I)=1
				ENDIF
			ELSEIF (P2=='F') THEN
				READ (IANS,'(X,I7,X,A2,3X,F16.9)') I,P1,FORCE
				IF (FORCE/=0) NUMC=NUMC+1
				IF (P1=='FX') THEN
					F(1,I)=FORCE
				ELSEIF (P1=='FY') THEN
					F(2,I)=FORCE
				ELSEIF (P1=='FZ') THEN
					F(3,I)=FORCE
				ENDIF
	        ELSE
				GOTO 1000
			ENDIF
		ENDDO
	ELSE
	    CYCLE
	ENDIF
ENDDO
!人工输入材料性质
1000WRITE(*,*) '输入单元类型 NPAR1=:(truss=1,4Q=2,3T=3,8H=4,beam=5,plate=6,shell=7) '
READ(*,*) NPAR1
WRITE(*,*) '输入杨氏模量 E=: '
READ(*,*) MAT%E
WRITE(*,*) '输入密度 DENS=: '
READ(*,*) MAT%DENS
IF (ELETY==180) THEN
    WRITE(*,*) '输入截面积 AREA=: '
    READ(*,*) MAT%AREA
ELSEIF (ELETY==182) THEN
    WRITE(*,*) '输入泊松比 PO=: '
    READ(*,*) MAT%P
ENDIF
!输出STAP90.IN
!输出标题行
WRITE (IIN,'(A40)') 'CESHI FOR YOUXIANYUAN---BY ZLX'     
!输出控制行        
WRITE (IIN,'(5I5)') NUMP,1,1,1,3       
!输出节点数据
DO I=1,NUMP
	WRITE (IIN,'(4I5,3F10.5)') I,ID(1,I),ID(2,I),ID(3,I),POS(1,I),POS(2,I),POS(3,I)
ENDDO
!输出载荷数据控制行
WRITE (IIN,'(2I5)') 1,NUMC
!输出各工况数据
DO I=1,NUMP
	DO J=1,3
		IF (F(J,I)/=0) THEN
			WRITE (IIN,'(2I5,F10.0)') I,J,F(J,I)
		ENDIF
	ENDDO
ENDDO
!输出单元控制行
IF (ELETY==180) THEN
    WRITE (IIN,'(4I5)') 1,NUME,1,1
ELSEIF (ELETY==182) THEN
    IF (NPAR1==2) THEN
        WRITE (IIN,'(4I5)') 2,NUME,1,1
    ELSE
        WRITE (IIN,'(4I5)') 3,2*NUME,1,1
    ENDIF
ENDIF

!输出材料数据
IF (ELETY==180) THEN
	WRITE (IIN,'(I5,E10.3,2F10.5)') 1,MAT%E,MAT%AREA,MAT%DENS
ELSEIF (ELETY==182) THEN
    WRITE (IIN,'(I5,E10.3,2F10.5)') 1,MAT%E,MAT%P,MAT%DENS
ENDIF
!输出单元数据
DO I=1,NUME
    IF (ELETY==180) THEN
	    WRITE (IIN,'(4I5)') I,IEN(1,I),IEN(2,I),1
    ELSEIF (ELETY==182) THEN
        IF (NPAR1==2) THEN
	        WRITE (IIN,'(6I5)') I,IEN(1,I),IEN(2,I),IEN(3,I),IEN(4,I),1
        ELSE
            WRITE (IIN,'(5I5)') 2*I-1,IEN(1,I),IEN(2,I),IEN(3,I),1
            WRITE (IIN,'(5I5)') 2*I,IEN(3,I),IEN(4,I),IEN(1,I),1
        ENDIF
	ENDIF
ENDDO

    END SUBROUTINE
    
