! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz

  
  !FTensDim = 0
  !RETURN
  
  f = 0
  g = 0
  h = 0
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE Parameters, ONLY :  nVar, nDim
   IMPLICIT NONE
   REAL, PARAMETER :: epsilon = 1e-14
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, nDim)
   REAL, INTENT(IN)  :: Q(nVar)
  ! Linear elasticity variables
   REAL :: lam,mu,irho, ialpha, u(3)
   REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar) 
   
  ! BgradQ = 0
  ! RETURN
  Qx = gradQ(:,1)
  Qy = gradQ(:,2)
  IF(nDim==3) THEN
	Qz = gradQ(:,3)
  ELSE
	Qz = 0.0 
  ENDIF 
  !
    ! Linear elasticity part
   AQx = 0.   
   BQy = 0. 
   CQz = 0. 

    lam  = Q(10)
    mu   = Q(11) 
    irho = 1./Q(12)
    IF(Q(13)<=1e-3) THEN
        BgradQ = 0.0  
        RETURN 
        !
    ELSE
        ialpha = 1./Q(13)
        u      = Q(7:9)*ialpha 
    ENDIF 
    !
    AQx(1) = - (lam+2*mu)*Qx(7) + (lam+2*mu)*u(1)*Qx(13) 
    AQx(2) = - lam*Qx(7)        + lam*u(1)*Qx(13)
    AQx(3) = - lam*Qx(7)        + lam*u(1)*Qx(13) 
    AQx(4) = - mu *Qx(8)        + mu *u(2)*Qx(13) 
    AQx(5) =   0.0 
    AQx(6) = - mu *Qx(9)        + mu *u(3)*Qx(13) 
    AQx(7) = - irho * Qx(1) - 2*Q(1)*irho*Qx(13)  
    AQx(8) = - irho * Qx(4) - 2*Q(4)*irho*Qx(13)   
    AQx(9) = - irho * Qx(6) - 2*Q(6)*irho*Qx(13)       
    !
    BQy(1) = - lam*Qy(8)        + lam*u(2)*Qy(13) 
    BQy(2) = - (lam+2*mu)*Qy(8) + (lam+2*mu)*u(2)*Qy(13) 
    BQy(3) = - lam*Qy(8)        + lam*u(2)*Qy(13)  
    BQy(4) = - mu *Qy(7)        + mu *u(1)*Qy(13) 
    BQy(5) = - mu *Qy(9)        + mu *u(3)*Qy(13) 
    BQy(6) =   0.0  
    BQy(7) = - irho * Qy(4) - 2*Q(4)*irho*Qy(13)  
    BQy(8) = - irho * Qy(2) - 2*Q(2)*irho*Qy(13)   
    BQy(9) = - irho * Qy(5) - 2*Q(5)*irho*Qy(13)      
    !
    CQz(1) = - lam*Qz(9)        + lam*u(3)*Qz(13) 
    CQz(2) = - lam*Qz(9)        + lam*u(3)*Qz(13) 
    CQz(3) = - (lam+2*mu)*Qz(9) + (lam+2*mu)*u(3)*Qz(13) 
    CQz(4) =  0.0  
    CQz(5) = - mu *Qz(8)        + mu *u(2)*Qz(13)  
    CQz(6) = - mu *Qz(7)        + mu *u(1)*Qz(13) 
    CQz(7) = - irho * Qz(6) - 2*Q(6)*irho*Qz(13)  
    CQz(8) = - irho * Qz(5) - 2*Q(5)*irho*Qz(13)  
    CQz(9) = - irho * Qz(3) - 2*Q(3)*irho*Qz(13)     

    
    BgradQ = AQx + BQy + CQz
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE Parameters, ONLY :  nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), V(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  REAL  :: lam,mu,irho
  ! Local Variables 
  REAL :: Pi
  Pi = ACOS(-1.0)

  L(:) = 0
  
    lam  = Q(10)   
    mu   = Q(11) 
    irho = 1./Q(12)
    !
    L(1) = -SQRT((lam+2.0*mu)*irho)     
    L(2) = +SQRT((lam+2.0*mu)*irho) 
    L(3) = -SQRT(mu*irho) 
    L(4) = +SQRT(mu*irho) 
    L(5) = 0. 

END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE Parameters, ONLY:  nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
  ! Local variable declaration 

  
  S = 0
      
END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(Name) 
  USE Parameters, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: Name(nVar)

  ! EQNTYPE99
    Name(1)  = 'sxx' 
    Name(2)  = 'syy' 
    Name(3)  = 'szz'
    Name(4)  = 'sxy' 
    Name(5)  = 'syz' 
    Name(6)  = 'sxz' 
    Name(7)  = 'u' 
    Name(8)  = 'v' 
    Name(9)  = 'w' 
    Name(10) = 'lambda'
    Name(11) = 'mu'
    Name(12) = 'rho' 
    Name(13) = 'alpha'
    Name(14) = 'xi' 
END SUBROUTINE PDEVarName


RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
  USE Parameters, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(nDim) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
    ! Linear elasticity variables
   REAL :: lam,mu,irho, ialpha, uv(3)
   REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar)
   !An = 0
  !RETURN

    lam  = Q(10) 
    mu   = Q(11) 
    irho = 1./Q(12)
    
    IF(Q(13)<=1e-3) THEN        
        ialpha = 0.0
        uv = 0.0 
        An = 0.0
        RETURN 
    ELSE
        ialpha = 1./Q(13)
        uv     = Q(7:9)*ialpha 
    ENDIF 
    
    A = 0.0
    B = 0.0
    C = 0.0 

    A(1,7)  = - (lam+2*mu) 
    A(1,13) = + (lam+2*mu)*uv(1)  
    A(2,7)  = - lam 
    A(2,13) = + lam*uv(1) 
    A(3,7)  = - lam 
    A(3,13) = + lam*uv(1)  
    A(4,8)  = - mu
    A(4,13) = + mu*uv(2)
    A(6,9)  = - mu
    A(6,13) = + mu*uv(3)
    A(7,1)  = - irho
    A(7,13) = - 2*Q(1)*irho  
    A(8,4)  = - irho
    A(8,13) = - 2*Q(4)*irho    
    A(9,6)  = - irho 
    A(9,13) = - 2*Q(6)*irho        
    !
    B(1,8)  = - lam 
    B(1,13) = + lam*uv(2) 
    B(2,8)  = - (lam+2*mu) 
    B(2,13) = + (lam+2*mu)*uv(2) 
    B(3,8)  = - lam 
    B(3,13) = + lam*uv(2) 
    B(4,7)  = - mu  
    B(4,13) = + mu*uv(1)  
    B(5,9)  = - mu 
    B(5,13) = + mu*uv(3) 
    B(7,4)  = - irho 
    B(7,13) = - 2*Q(4)*irho 
    B(8,2)  = - irho 
    B(8,13) = - 2*Q(2)*irho 
    B(9,5)  = - irho 
    B(9,13) = - 2*Q(5)*irho       
    !
    C(1,9)  = - lam 
    C(1,13) = + lam*uv(3) 
    C(2,9)  = - lam 
    C(2,13) = + lam*uv(3) 
    C(3,9)  = - (lam+2*mu) 
    C(3,13) = + (lam+2*mu)*uv(3)
    C(5,8)  = - mu  
    C(5,13) = + mu*uv(2)  
    C(6,7)  = - mu 
    C(6,13) = + mu*uv(1) 
    C(7,6)  = - irho 
    C(7,13) = - 2*Q(6)*irho 
    C(8,5)  = - irho 
    C(8,13) = - 2*Q(5)*irho 
    C(9,3)  = - irho 
    C(9,13) = - 2*Q(3)*irho     
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if
    
    

  
END SUBROUTINE PDEMatrixB


