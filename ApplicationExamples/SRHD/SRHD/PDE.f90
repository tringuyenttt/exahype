! Special Relativistic Hydrodynamics Conservative Formulation

SUBROUTINE PDEFlux(F,Q) 
  USE Parameters, ONLY : nVar, nDim, gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  REAL, PARAMETER :: epsilon = 1e-14
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar) 
  REAL, INTENT(OUT) :: F(nVar,nDim) 
  ! Local variables  
  REAL :: rho,vx,vy,vz,p
  REAL :: v2,lf,w,ww,wwx,wwy,wwz,gamma1
  REAL :: V(nVar)
  INTEGER :: iErr

  CALL PDECons2Prim(V,Q,iErr)
  gamma1 = gamma/(gamma-1.0)
  rho    = V(1)
  vx     = V(2)
  vy     = V(3)
  vz     = V(4)
  p      = V(5)

  v2     = vx**2 + vy**2 + vz**2
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p
  ww     = w*lf**2
  wwx    = ww*vx
  wwy    = ww*vy
  wwz    = ww*vz

  F(1, 1) = vx*rho*lf
  F(2, 1)   = wwx*vx + p
  F(3, 1)   = wwx*vy 
  F(4, 1)   = wwx*vz  
  F(5, 1)   = wwx - F(1, 1)

  F(1, 2)   = vy*rho*lf
  F(2, 2)   = wwy*vx  
  F(3, 2)   = wwy*vy + p
  F(4, 2)   = wwy*vz  
  F(5, 2)   = wwy - F(1, 2)

  IF ( nDim > 2 ) then
    F(1, 3)   = vz*rho*lf
    F(2, 3)   = wwz*vx  
    F(3, 3)   = wwz*vy  
    F(4, 3)   = wwz*vz + p
    F(5, 3)   = wwz - F(1, 3)
  END IF

END SUBROUTINE PDEFlux 



SUBROUTINE PDEEigenvalues(Lambda,Q,nv) 
  USE Parameters, ONLY : nVar, nDim, gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar), nv(nDim)
  REAL, INTENT(OUT) :: Lambda(nVar)  
  ! Local variables  
  INTEGER :: iErr
  REAL :: rho, vel(nDim), p
  REAL :: cs2, cs, c0, v2, w, gamma1, vn, den, u, c
  REAL :: V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14  
  Lambda(:) = 0.0
  !
  CALL PDECons2Prim(V,Q,iErr)
  rho    = V(1)
  ! the velocity vector is dimension agnostic
  vel    = V(2:2+nDim-1)
  p      = V(5)
  gamma1 = gamma/(gamma-1.0)
  w      = rho + gamma1*p
  cs2    = gamma*p/w
  v2     = SUM(vel**2)
  vn     = SUM(vel * nv)
  den    = 1.0/(1.0 - v2*cs2)
  IF(SUM(nv**2).EQ.0.) THEN  
     u = SQRT( v2) 
  ELSE
     u = vn 
  ENDIF
  Lambda(1)   = ( u*(1.0-cs2)-SQRT( cs2*(1.0-v2)*( (1.0-v2*cs2) - u**2*(1.0-cs2) )) )*den
  Lambda(2:4) = u
  Lambda(5)   = ( u*(1.0-cs2)+SQRT( cs2*(1.0-v2)*( (1.0-v2*cs2) - u**2*(1.0-cs2) )) )*den 
  !
END SUBROUTINE PDEEigenvalues
 

