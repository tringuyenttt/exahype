!!!
!!! A Special Relativistic ideal Magnetohydrodynamics ExaHyPE Kernel.
!!! Based on files provided by Olindo, adopted by Sven.
!!! Based on srhd3dfortran/PDE.f90.
!!!

SUBROUTINE PDEEigenvalues(Lambda,Q,nv) 
  USE Parameters, ONLY : nVar, nDim , gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar), nv(nDim)  
  REAL, INTENT(OUT) :: Lambda(nVar)  
  ! Local variables  
  INTEGER :: iErr
  REAL :: rho,vx,vy,vz,p
  REAL :: cs2, cs, c0, v2, w, gamma1, vn, den, u, c
  REAL :: ex, ey, ez, b2, e2, a2, ca2, vf1, vf2, vel(nDim)
  REAL :: bx, by, bz
  REAL :: V(nVar)

  ! These are not the exact eigenvalues, instead of Lambda(1..9)
  ! we compute only two eigenvalues: Approximate magnetosonics

  CALL PDECons2Prim(V,Q,iErr)
  rho    = V(1)
  vx     = V(2)
  vy     = V(3)
  vz     = V(4)
  p      = V(5)
  bx     = V(6)
  by     = V(7)
  bz     = V(8)

  ex     = - (vy*bz - vz*by)
  ey     = - (vz*bx - vx*bz)
  ez     = - (vx*by - vy*bx)
  gamma1 = gamma/(gamma-1.0)
  b2     = bx*bx + by*by + bz*bz
  !v2     = vx*vx + vy*vy + vz*vz
  e2     = ex*ex + ey*ey + ez*ez
  w      = rho + gamma1*p
  !vn     = vx*nv(1) + vy*nv(2) + vz*nv(3)

  ! the velocity vector is dimension agnostic
  vel    = V(2:2+nDim-1)
  v2     = SUM(vel**2)
  vn     = SUM(vel * nv)

  cs2    = gamma * p / w
  ca2    = (b2 - e2) / ( w + b2 - e2)
  a2     = cs2 + ca2 - cs2 * ca2
  den    = 1.0/( 1.0 - v2 * a2)
  vf1    = den * ( 1.0 - a2) * vn
  vf2    = den * sqrt(a2 * ( 1.0 - v2) * ((1.0 - v2 * a2) - (1.0 - a2) * vn**2))

  Lambda(1) = vf1 + vf2
  Lambda(2) = vf1 - vf2
  Lambda(3) = 0.0
  Lambda(4) = 0.0
  Lambda(5) = 0.0
  Lambda(6) = 0.0
  Lambda(7) = 0.0
  Lambda(8) = 0.0
  Lambda(9) = 0.0

END SUBROUTINE PDEEigenvalues
 


SUBROUTINE PDEFlux(F,Q) 
  USE Parameters, ONLY : nVar, nDim, gamma, DivCleaning_a
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar) 
  REAL, INTENT(OUT) :: F(nVar,nDim) 
  ! Local variables  
  REAL :: rho,vx,vy,vz,p
  REAL :: v2,lf,w,ww,wwx,wwy,wwz,gamma1
  REAL :: ex,ey,ez,b2,e2,uem,bx,by,bz
  REAL :: V(nVar)
  INTEGER :: iErr
  
  CALL PDECons2Prim(V,Q,iErr)
  gamma1 = gamma/(gamma-1.0)
  rho    = V(1)
  vx     = V(2)
  vy     = V(3)
  vz     = V(4)
  p      = V(5)
  bx     = V(6)
  by     = V(7)
  bz     = V(8)

  ex     = - (vy*bz - vz*by)
  ey     = - (vz*bx - vx*bz)
  ez     = - (vx*by - vy*bx)

  v2     = vx**2 + vy**2 + vz**2
  b2     = bx**2 + by**2 + bz**2
  e2     = ex**2 + ey**2 + ez**2
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p
  ww     = w*lf**2
  uem    = 0.5*(b2 + e2)
  wwx    = ww*vx
  wwy    = ww*vy
  wwz    = ww*vz

  !
  F(1, 1) = vx*rho*lf
  F(2, 1) = wwx*vx - bx*bx - ex*ex + p + uem
  F(3, 1) = wwx*vy - bx*by - ex*ey
  F(4, 1) = wwx*vz - bx*bz - ex*ez 
  F(5, 1) = wwx + (ey*bz - ez*by) 
  F(6, 1) = V(9)
  F(7, 1) = -ez
  F(8, 1) = ey  
  F(9, 1) = DivCleaning_a**2*bx
  !
  F(1, 2) = vy*rho*lf
  F(2, 2) = wwy*vx - by*bx - ey*ex 
  F(3, 2) = wwy*vy - by*by - ey*ey + p + uem
  F(4, 2) = wwy*vz - by*bz - ey*ez 
  F(5, 2) = wwy + (ez*bx - ex*bz) 
  F(6, 2) = ez 
  F(7, 2) = V(9) 
  F(8, 2) = -ex   
  F(9, 2) = DivCleaning_a**2*by
  !
  IF ( nDim > 2 ) then
    F(1, 3) = vz*rho*lf
    F(2, 3) = wwz*vx - bz*bx - ez*ex 
    F(3, 3) = wwz*vy - bz*by - ez*ey 
    F(4, 3) = wwz*vz - bz*bz - ez*ez + p + uem
    F(5, 3) = wwz + (ex*by - ey*bx) 
    F(6, 3) = -ey  
    F(7, 3) = ex   
    F(8, 3) = V(9)   
    F(9, 3) = DivCleaning_a**2*bz
  ENDIF
  !
END SUBROUTINE PDEFlux 
 
