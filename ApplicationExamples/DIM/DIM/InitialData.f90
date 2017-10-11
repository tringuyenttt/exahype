! DIM Initial Data


RECURSIVE SUBROUTINE InitialData(x, t, Q)
	USE, INTRINSIC :: ISO_C_BINDING
	USE Parameters, ONLY : nVar, nDim
	IMPLICIT NONE 
	! Argument list 
	REAL, INTENT(IN)               :: x(nDim), t        ! 
	REAL, INTENT(OUT)              :: Q(nVar)        ! 

	Call InitialPlaneWave(x, t, Q)
    	!call GaussianBubble(x,t,Q)
	END SUBROUTINE InitialData


RECURSIVE SUBROUTINE InitialPlaneWave(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r
    ! Initialize parameters
    ICuL(:)=(/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0, 1.0 /)
    ICuR(:)=(/ 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
    ICA=1.0     ! wave length
    ICxd=0.25   ! radius of the circle
    ICsig=1e-6  ! smoothing parameter
    ICQR(:)= (/1e-7, 0.9999999 /)
    ! Initialize the variables vector V
    up = ICuL + ICuR*SIN(2*Pi*(x(1)-2.0*t)/ICA)
    r = SQRT(sum(x(1:nDim)**2)) 
    xi = 0.5+0.5*ERF((r-ICxd)/ICsig) 
    !up(13)  = ICQR(1)*(1-xi) + ICQR(2)*xi 
    up(13)=1.0
    up(1:9) = up(1:9)*up(13)
    up(14)=1.0-1.e-10
    
    !CALL PDEPrim2Cons(Q,V)
    Q=up
    END SUBROUTINE InitialPlaneWave
    
RECURSIVE SUBROUTINE GaussianBubble(x, t, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL    :: up(nVar), Pi = ACOS(-1.0)
    REAL    :: xi, ICQR(2), ICxd, ICsig
    REAL    :: ICuL(nVar), ICuR(nVar),ICA,r

        up=0;
        r = SQRT((x(1))**2 + (x(2))**2)
        up(1:2)=exp(-40.0*r**2)
        up(13)=1.0
        up(14)=1.0
        up(10) = 2.0  
        up(11) = 1.0  
        up(12) = 1.0 
    
    !CALL PDEPrim2Cons(Q,up)
        Q=up
END SUBROUTINE GaussianBubble

