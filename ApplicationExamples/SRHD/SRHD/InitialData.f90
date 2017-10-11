 ! 
 ! This file is part of the ExaHyPE project.
 ! Copyright (c) 2016  http://exahype.eu
 ! All rights reserved.
 !
 ! The project has received funding from the European Union's Horizon 
 ! 2020 research and innovation programme under grant agreement
 ! No 671698. For copyrights and licensing, please consult the webpage.
 !
 ! Released under the BSD 3 Open Source License.
 ! For the full license text, see LICENSE.txt
 !  
 
 SUBROUTINE MinimumTreeDepth(depth)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    INTEGER, INTENT(OUT)              :: depth        ! maximal depth of tree recursion
    depth = 4
END SUBROUTINE MinimumTreeDepth 

SUBROUTINE AdjustedSolutionValues(x, w, t, dt, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(2)        ! 
	REAL, INTENT(IN)               :: w           ! 
	REAL, INTENT(IN)               :: t           ! 
	REAL, INTENT(IN)               :: dt          ! 

	REAL, INTENT(OUT)              :: Q(5)        ! 
	
	IF ( t < 1e-15 ) THEN
		CALL InitialData(x, Q)
	ENDIF
END SUBROUTINE AdjustedSolutionValues


SUBROUTINE InitialData(x, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(2)        ! 
	REAL, INTENT(OUT)              :: Q(5)        ! 

	! this implements the 1D Riemann problem as initial data

	REAL, PARAMETER :: gamma=5./3., t_f=0.4, xsep=0.5
	REAL, PARAMETER :: rho_L= 1, v_L=-0.6, p_L=10
	REAL, PARAMETER :: rho_R=10, v_R=+0.5, p_R=20
	REAL :: V(nVar)

	IF (x(1) < xsep) THEN
	        V(1) = rho_L
	        V(2) = v_L
	        V(3) = 0.
	        V(4) = 0.
	        V(5) = p_L
	ELSE
	        V(1) = rho_R
	        V(2) = v_R
	        V(3) = 0.
	        V(4) = 0.
	        V(5) = p_R
	END IF

	Call PDEPrim2Cons(Q,V)
END SUBROUTINE InitialData




