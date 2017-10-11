! The Con2Prim and Prim2Con routines for SRHD.
! Should be merged with MHD's.

SUBROUTINE PDECons2Prim(V,Q,iErr)
  USE parameters, ONLY : nVar, gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  !REAL :: x(3), time
 ! Removed as not present and also not used
  INTEGER :: iErr
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variable declaration
  REAL        :: iRho
  REAL, PARAMETER :: epsilon = 1e-14
  INTEGER     :: i, iter, indx_nr(2)
  REAL        :: sx, sy, sz, bx, by, bz
  REAL        :: v2, sb, den, vb, zeta, cs    
  REAL        :: gamma1, G1, G12, x1, x2, eps
  REAL        :: rho, vx, vy, vz, p
  REAL        :: lf, lf2, lf3, lf4, gam,d,e,s2,b2,e2,sb2,w,ww
  REAL        :: RTSAFE_C2P_RMHD1, RTSAFE_C2P_RMHD2, RTSAFE_C2P_RHD1, RTSAFE_C2P_RHD2
  REAL        :: k(3), B(3), vel(3)
  REAL        :: p1,q1,dd,phi,temp1, H, k2, kB, T2
  REAL        :: epsilon0, mu0, f, df, drho
  REAL        :: iPi,B28P, dcs, dc0
  REAL        :: dv2, c0, dw
  REAL        :: ri,qi,kappa,z,kappa_max
  REAL        :: ZBRENT_C2P_RHD2, ZBRENT_LF_POLY
  REAL        :: lapse, gp, gm
  REAL        :: g_contr(3,3), g_cov(3,3)
  REAL        :: shift(3), sm(3), sm_cov(3), vf_cov(3), vf(3), Ev(3), Bv(3)
  REAL        :: Qtest(14), dQ(14)
  REAL        :: ExB(3),ExB_up(3), U2e(3), S_up(3)  
  REAL        :: LFsolutions(4), cLF(5), C2, C3 
  COMPLEX     :: LFsolutionsC(4)
  REAL        :: xi(2), alpha_nr(2,2), beta_nr(2), d_nr, S  
  REAL        :: A(3,3), devG(3,3), G(3,3), tempp(3,3), Id(3,3), detA, ehh, evv
  REAL        :: EE(3), vv(3), BB(3), eel, detvEB, vxE(3), vxB(3)  
  REAL        :: vs,vl,us,ul,rhos,rhol,rhoeh,rhoe,alphas,alphal
  REAL        :: xx(9), temp(3) 
  LOGICAL     :: FAILED
  REAL, PARAMETER    :: tol = 1e-8, third=1.0/3.0, p_floor = 1.0e-5, rho_floor = 1.0e-4

  iErr = 0

  gamma1 = gamma/(gamma - 1.0)
  gam    = 1.0/gamma1

  IF (Q(1) .LT. 0.0) THEN
     rho = rho_floor
     vx  = 0.0
     vy  = 0.0
     vz  = 0.0
     p   = p_floor
     V(1:nVar) = (/ rho, vx, vy, vz, p /)
     RETURN
  ENDIF
  ! First option: solve for v2 (Del Zanna et al. (2007) A&A, 473, 11-30)
  FAILED = .FALSE.
  d    = Q(1)
  sx   = Q(2)
  sy   = Q(3)
  sz   = Q(4)
  e    = Q(5) + d
  bx = 0.0
  by = 0.0
  bz = 0.0
  !
  s2   = sx*sx + sy*sy + sz*sz
  b2   = 0.0
  sb   = 0.0
  sb2  = 0.0
  eps  = 1.e-8
  !
  x1   = 0.
  x2   = 1.-eps
  v2   = RTSAFE_C2P_RMHD1(x1,x2,tol,gam,d,e,s2,b2,sb2,w,FAILED)

  IF (FAILED) THEN
     iErr = -1
     p    = p_floor
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
  ELSE

     den  = 1.0/(w+b2)

     vb   = sb/w
     !
     rho  = d*sqrt(1.-v2)
     vx   = sx*den
     vy   = sy*den
     vz   = sz*den
     p    = max(1.e-15, gam*(w*(1.-v2)-rho))
     !
  ENDIF
  !! Second option: solve for p [Appendix D of the book "Relativistic hydrodynamics" by Rezzolla & Zanotti (2013)]
  !FAILED = .FALSE.
  !d    = Q(1)
  !sx   = Q(2)
  !sy   = Q(3)
  !sz   = Q(4)
  !e    = Q(5) 
  !!
  !s2   = sx*sx + sy*sy + sz*sz
  !eps  = 1.e-8
  !!
  !x1   = eps      ! min pressure
  !x2   = 1.0d+5   ! max pressure
  !p    = RTSAFE_C2P_RHD1(x1,x2,tol,d,e,s2,FAILED)
  !
  !IF (FAILED) THEN
  !   iErr = -1
  !   p    = p_floor
  !   rho  = rho_floor
  !   vx   = 0.0
  !   vy   = 0.0
  !   vz   = 0.0
  !ELSE
  !   rho  = d / (e + p + d) * SQRT((e + p + d)**2 - s2)
  !   den  = 1.0 / (e + p + d)
  !   !
  !   vx   = sx*den
  !   vy   = sy*den
  !   vz   = sz*den
  !ENDIF

  ! Third option: solve for z = lf*v [Galeazzi et al. (2013), Phys Rev D 88 064009, Appendix C]
  !
  !  FAILED = .FALSE.
  !  d    = Q(1)
  !  sx   = Q(2)
  !  sy   = Q(3)
  !  sz   = Q(4)
  !  e    = Q(5) 
  !  !
  !  s2   = sx*sx + sy*sy + sz*sz
  !  ri   = SQRT(s2)/d
  !  qi   = e/d
  !  kappa= ri/(1.0 + qi)
  !  kappa_max = 0.999  
  !  IF (kappa > kappa_max) THEN
  !    ri   = kappa_max*(1.0 + qi)
  !    kappa= kappa_max
  !  ENDIF
  !  x1   = 0.5*kappa/SQRT(1.0-0.25*kappa**2)    ! min z
  !  x2   =     kappa/SQRT(1.0-     kappa**2)    ! max z
  !!  z    = ZBRENT_C2P_RHD2 ( x1, x2, tol, qi, ri, kappa, FAILED)
  !  z    = RTSAFE_C2P_RHD2(x1, x2, tol, qi, ri, kappa, FAILED)
  !  !
  !  IF (FAILED) THEN
  !    iErr = -1
  !    p    = p_floor
  !    rho  = rho_floor
  !    vx   = 0.0
  !    vy   = 0.0
  !    vz   = 0.0
  !  ELSE
  !    lf    = SQRT(1.0 + z**2)
  !    rho   = d/lf
  !    eps   = lf*(qi+1.0) - z*ri - 1.0
  !
  !    w     = (1.0 + gamma*eps)*d*lf
  !    den   = 1.0/w
  !  
  !    p     = rho*eps*(gamma - 1.0)
  !    vx    = sx*den
  !    vy    = sy*den
  !    vz    = sz*den
  !  ENDIF
  !
  ! Fourth option: Isentropic case, for which the energy equation is not needed.
  !FAILED = .FALSE.
  !d    = Q(1)
  !sx   = Q(2)
  !sy   = Q(3)
  !sz   = Q(4)
  !
  !s2 = sx*sx + sy*sy + sz*sz
  !b2 = 0.0
  !sb = 0.0
  !sb2 = sb**2
  !
  !x1  = 1.0
  !x2  = 10.0
  !
  !! Compute the Lorentz factor 
  !lf = ZBRENT_LF_POLY ( x1, x2, tol, gamma, IC%IsenRHD%KK, d, b2, s2, sb)
  !
  !rho = d / lf
  !IF ( rho < rho_floor ) THEN
  !   rho = rho_floor
  !ENDIF
  !
  !! The following expression is for an isentropic polytrope only
  !p     = IC%IsenRHD%KK * rho**gamma
  !eps   = IC%IsenRHD%KK / (gamma - 1.0) * rho**(gamma - 1.0)
  !
  !h     = 1.0d0 + gamma * eps
  !zeta  = rho * h * lf * lf
  !den   = 1.0 / zeta
  !
  !vx = sx * den
  !vy = sy * den
  !vz = sz * den  
  !
  V(1:nVar) = (/ rho, vx, vy, vz, p /)
END SUBROUTINE PDECons2Prim


! This subroutine is not used in this fortran code but only by the initial conditions,
! which are given in C++ code. We need an interface for that.
SUBROUTINE PDEPrim2Cons(Q,V)
  USE parameters, ONLY : nVar, gamma
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
  REAL :: Prim(nVar), Buf(nVar)
  REAL :: rho,vx,vy,vz,p,bx,by,bz,ex,ey,ez,cs,c0
  REAL :: v2,b2,e2,lf,w,ww,uem,gamma1
  REAL :: lapse, gp, gm, dcs, dc0, eel 
  REAL :: g_contr(3,3), g_cov(3,3)
  REAl :: shift(3), vf(3), vf_cov(3), Ev(3), Bv(3), ExB(3)
  REAL :: A(3,3), devG(3,3), G(3,3), temp(3,3), Id(3,3), detA, eh, S, evv, T, falpha  
  REAL :: alphas, alphal, rhos, rhol, us, vs, ul, vl, es, el
  REAL :: EE(3),BB(3),vv(3),detvEB, vxE(3), vxB(3)  
  INTEGER :: i

  !
  rho    = V(1)
  vx     = V(2)
  vy     = V(3)
  vz     = V(4)
  p      = V(5)
  !
  v2     = vx**2 + vy**2 + vz**2
  !
  IF (v2 > 1.0) THEN
     WRITE(*,*)'Superluminal velocity in PDEPrim2Cons!!'
     STOP
  ENDIF
  lf     = 1.0 / sqrt(1.0 - v2)
  gamma1 = gamma/(gamma-1.0)
  w      = rho + gamma1*p
  ww     = w*lf**2
  !
  Q(1)   = rho*lf
  Q(2)   = ww*vx 
  Q(3)   = ww*vy 
  Q(4)   = ww*vz 
  Q(5)   = ww - p - Q(1)
END SUBROUTINE PDEPrim2Cons
