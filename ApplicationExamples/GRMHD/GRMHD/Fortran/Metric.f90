! Static background metrics for GRMHD

RECURSIVE SUBROUTINE METRIC ( xc, lapse, gp, gm, shift, g_cov, g_contr)
  USE Parameters, ONLY : nVar, nDim
  IMPLICIT NONE
  !
  REAL, DIMENSION(nDim), intent(IN) :: xc
  REAL                :: lapse,gp,gm
  REAL,dimension(3)   :: shift
  REAL,dimension(3,3) :: g_cov, g_contr

  REAL :: x, y, z, r, z2, r2, aom2
  REAL :: lx, ly, lz, HH, SS, detg
  REAL :: st, st2, delta, rho2, sigma, zz
  
  REAL :: aom = 0.0
  REAL :: Mbh = 1.0

  !

  ! Rotating black hole in Kerr-Schild spherical coordinates. See Appendix B of Komissarov (2004) MNRAS, 350, 427
  r  = xc(1)
  r2 = r*r
 
  st = SIN(xc(2))
  IF (st < 1.e-6) st = 1.
  st2 = st**2
 
  aom2 = aom**2
 
  delta = r2 - 2.0 * Mbh * r + aom2
  rho2  = r2 + aom2*(1.0 - st2)
  sigma = (r2 + aom2)**2 - aom2*delta*st2
  zz    = 2.0*r/rho2
 
  lapse    = 1.0 / sqrt(1.0 + zz)
  shift(1) = zz/(1.0 + zz)      
  shift(2) = 0.0        
  shift(3) = 0.0 
 
  g_cov( 1, 1:3) = (/ 1.0 + zz,             0.0,        -aom*st2*(1.0 + zz)   /)
  g_cov( 2, 1:3) = (/ 0.0,                  rho2,       0.0                   /)
  g_cov( 3, 1:3) = (/ -aom*st2*(1.0 + zz),  0.0,        (sigma/rho2)*st2      /)
 
  g_contr( 1, 1:3) = (/  lapse**2 + aom2*st2/rho2,   0.0,                       aom/rho2         /)
  g_contr( 2, 1:3) = (/  0.0,                        1.0/rho2,                  0.0              /)
  g_contr( 3, 1:3) = (/  aom/rho2,                   0.0,                       1.0/(rho2*st2)   /)
 
  gp = rho2*st*sqrt(1.0 + zz)
  gm = 1.0/gp 
  !
 
END SUBROUTINE METRIC

SUBROUTINE METRIC_3D  ( xc, lapse, gp, gm, shift, g_cov, g_contr)
  USE Parameters, ONLY : nVar, nDim
  IMPLICIT NONE
  !
  REAL, DIMENSION(nDim), intent(IN) :: xc
  REAL                :: lapse,gp,gm
  REAL,dimension(3)   :: shift
  REAL,dimension(3,3) :: g_cov, g_contr

  REAL :: x, y, z, r, z2, r2, aom2
  REAL :: lx, ly, lz, HH, SS, detg
  REAL :: st, st2, delta, rho2, sigma, zz
  
  REAL :: aom = 0.0
  REAL :: Mbh = 1.0

! Rotating black hole in Kerr-Schild Cartesian coordinates. See De Felice & Clarke Sect. 11.4
  x    = xc(1)
  y    = xc(2)
  z    = xc(3)
 
  z2   = z**2
  aom2 = aom**2
  r    = SQRT( (x**2 + y**2 + z**2 - aom2)/2.0 + SQRT(((x**2 + y**2 + z**2 - aom2)/2.0)**2 + z2*aom2))
 
  r2   = r**2
 
  HH = Mbh*r2*r / (r2*r2 + aom2*z2)
  SS = 1.0 + 2.0*HH
  lx = (r*x + aom*y)/(r2 + aom2) 
  ly = (r*y - aom*x)/(r2 + aom2)
  lz = z/r 
 
  lapse   = 1.0/SQRT(SS)
  shift(1) = 2.0*HH/SS*lx        
  shift(2) = 2.0*HH/SS*ly        
  shift(3) = 2.0*HH/SS*lz       
    
  g_cov( 1, 1:3) = (/ 1.0 + 2.0*HH*lx**2, 2.0*HH*lx*ly,        2.0*HH*lx*lz       /)
  g_cov( 2, 1:3) = (/ 2.0*HH*lx*ly,       1.0 + 2.0*HH*ly**2,  2.0*HH*ly*lz       /)
  g_cov( 3, 1:3) = (/ 2.0*HH*lx*lz,       2.0*HH*ly*lz,        1.0 + 2.0*HH*lz**2 /)
 
  g_contr( 1, 1:3) = (/  1.0 + 2.0*HH*ly**2 + 2.0*HH*lz**2,   -2.0*HH*lx*ly,    -2.0*HH*lx*lz   /)
  g_contr( 2, 1:3) = (/ -2.0*HH*lx*ly,        1.0 + 2.0*HH*lx**2 + 2.0*HH*lz**2,  -2.0*HH*ly*lz    /)
  g_contr( 3, 1:3) = (/ -2.0*HH*lx*lz,     -2.0*HH*ly*lz,       1.0 + 2.0*HH*lx**2 + 2.0*HH*ly**2 /)
 
  g_contr = g_contr/SS
 
  gp = SQRT(SS)
  gm = 1.0/gp
END SUBROUTINE METRIC_3D
