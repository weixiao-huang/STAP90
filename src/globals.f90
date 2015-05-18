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

! .  Define global variables

module GLOBALS

   integer, parameter :: IELMNT=1	! Unit storing element data 
   integer, parameter :: ILOAD=2	! Unit storing load vectors
   integer, parameter :: IIN=5		! Unit used for input
   integer, parameter :: IOUT=6		! Unit used for output

   integer :: NUMNP		! Total number of nodal points
						! = 0 : Program stop
   integer :: NEQ		! Number of equations
   integer :: NWK		! Number of matrix elements
   integer :: MK		! Maximum half bandwidth

   integer :: IND		! Solution phase indicator
						!   1 - Read and generate element information
						!   2 - Assemble structure stiffness matrix
						!   3 - Stress calculations
   integer :: NPAR(10)	! Element group control data
						!   NPAR(1) - Element type
						!             1 : Truss element
						!   NPAR(2) - Number of elements
						!   NPAR(3) - Number of different sets of material and 
						!             cross-sectional  constants
                        !   NPAR(4) - 0-轴对称问题  1-平面应变  2-平面应力
                        !   NPAR(5) - 1-子空间迭代法  2-Lanczos法
   integer :: NUMEG		! Total number of element groups, > 0

   integer :: MODEX		! Solution mode: 0 - data check only;  1 -  execution  
   integer :: NDF       ! the degree of the freedom of the point in the element

   real :: TIM(5)		! Timing information
   character*80 :: HED	! Master heading information for use in labeling the output

   integer :: NFIRST
   integer :: NLAST
   integer :: NUMEST
   integer :: MIDEST
   integer :: MAXEST

   integer :: NG

end module GLOBALS
