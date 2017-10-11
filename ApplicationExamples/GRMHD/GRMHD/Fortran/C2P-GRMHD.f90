! C2P for GRMHD

RECURSIVE SUBROUTINE PDEPrim2Cons(Q,V)
  USE Parameters, ONLY: gamma, nVar, nDim 
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
  REAL :: psi, BV_contr(3), Qv_contr(3), QB_contr(3), vxB_contr(3), vb_conv(3), b2_cov, vb_cov
  INTEGER :: i


  rho     = V(1)
  vf_cov  = V(2:4)
  p       = V(5)
  !
  BV(1:3) = V(6:8)
  psi = V(9)
  lapse = V(10)
  shift = V(11:13)          ! NB: we choose V() and Q() being the shift_controvariant!
  !
  !gammaij = V(14:19) 
  g_cov(1,1) = V(14)
  g_cov(1,2) = V(15)
  g_cov(1,3) = V(16)
  g_cov(2,2) = V(17)
  g_cov(2,3) = V(18)
  g_cov(3,3) = V(19)
  g_cov(2,1) = V(15)
  g_cov(3,1) = V(16)
  g_cov(3,2) = V(18)
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  gp = SQRT(gp)
  gm = 1./gp
  !  
  !CALL METRIC(x, lapse, gp, gm, shift, g_cov, g_contr)
  !
  vf      = MATMUL(g_contr,vf_cov)
  BV_contr = MATMUL(g_contr,BV(1:3))
  Qv_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  vxB_contr = MATMUL(g_contr,vxB(1:3))
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  vb_cov      = vf_cov(1)*BV(1) + vf_cov(2)*BV(2) + vf_cov(3)*BV(3) 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  !
  b2_cov      = SUM(BV(1:3)**2)
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
  Q(1)    = rho*lf
  Q(2:4)  = ww*vf_cov(1:3) + b2_cov*vf_cov(1:3) - vb_cov*BV(1:3)
  Q(5)    = ww - p + uem - Q(1)     !!!!! we subtract PDE(Q(1))!!!!
  Q(6:8) = V(6:8)
  Q(1:8)    = gp*Q(1:8)
  Q(9:)    = V(9:)
!
END SUBROUTINE PDEPrim2Cons


RECURSIVE SUBROUTINE PDECons2Prim(V,Q,iErr)
  USE Parameters, ONLY: gamma, nVar, nDim
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTEGER :: iErr
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variable declaration
  REAL        :: iRho
  REAL, PARAMETER :: epsilon = 1e-14
  INTEGER     :: i, iter, indx_nr(2)
  REAL        :: dr, sx, sy, sz, bx, by, bz
  REAL        :: v2, sb, den, vb, zeta, cs    
  REAL        :: gamma1, G1, G12, x1, x2, eps
  REAL        :: rho, vx, vy, vz, p
  REAL        :: lf, lf2, lf3, lf4, gam,e,s2,b2,e2,sb2,w,ww
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
  REAL :: psi, BV_contr(3), Qv_contr(3), QB_contr(3), vxB_contr(3), vb_conv(3), b2_cov, vb_cov, gammaij(6), Qloc(nVar), d

  !
  iErr = 0

  !  
  psi = Q(9)
  lapse = Q(10)
  shift = Q(11:13)
  !
  gammaij = Q(14:19) 
  g_cov(1,1) = Q(14)
  g_cov(1,2) = Q(15)
  g_cov(1,3) = Q(16)
  g_cov(2,2) = Q(17)
  g_cov(2,3) = Q(18)
  g_cov(3,3) = Q(19)
  g_cov(2,1) = Q(15)
  g_cov(3,1) = Q(16)
  g_cov(3,2) = Q(18)
  !
  CALL MatrixInverse3x3(g_cov,g_contr,gp) 
  gp = SQRT(gp)
  gm = 1./gp
  ! 
  !CALL METRIC(x, lapse, gp, gm, shift, g_cov, g_contr)  
  Qloc(:)  = gm*Q(:)
  !
  BV(1:3) = Qloc(6:8) ! magnetic field
  !
  gamma1 = gamma/(gamma - 1.0)
  gam    = 1.0/gamma1
  ! Solve for p
  FAILED  = .FALSE.
  d       = Qloc(1)
  sm_cov  = Qloc(2:4)
  !
  sm   = MATMUL (g_contr, sm_cov)
  BV_contr   = MATMUL (g_contr, BV)
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  b2   = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  sb   = sm_cov(1)*BV_contr(1) + sm_cov(2)*BV_contr(2) + sm_cov(3)*BV_contr(3)
  sb2  = sb**2
  eps  = 1.e-10 !8
  !!! RHD
  !!e       = Qloc(5) 
  !!x1   = eps      ! min pressure
  !!x2   = 1.0d+5   ! max pressure
  !!p    = RTSAFE_C2P_RHD1(x1,x2,tol,d,e,s2,FAILED)
  !!IF (FAILED) THEN
  !!   p = 1.0e-20
  !!ENDIF
  !!rho  = d / (e + p + d) * SQRT((e + p + d)**2 - s2)
  !!den  = 1.0 / (e + p + d)
  !!!
  !!vf_cov(1:3) = sm_cov(1:3)*den
  !!
  ! First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
  e    = Qloc(5) + d  ! Q(5) = gamma^1/2 ( U - D )
  x1   = 0.      ! 
  x2   = 1.0-eps ! 
  w=0
  v2   = RTSAFE_C2P_RMHD1(x1,x2,tol,gam,d,e,s2,b2,sb2,w,FAILED)
  
  !PRINT *, "FC2P: v2=",v2,"w=",w
  
  !
  IF (FAILED) THEN
     iErr = -1
     p    = p_floor
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
     bx   = bx
     by   = by
     bz   = bz
  ELSE
     den  = 1.0/(w+b2)
     vb   = sb/w
     !
     rho  = d*sqrt(1.-v2)
     vf_cov(1) = (sm_cov(1) + vb*BV(1))*den
     vf_cov(2) = (sm_cov(2) + vb*BV(2))*den
     vf_cov(3) = (sm_cov(3) + vb*BV(3))*den
     p = max(1.e-15, gam*(w*(1.-v2)-rho))
  ENDIF
  !
  V(1:19) = (/ rho, vf_cov(1:3), p, BV(1:3), psi , lapse, shift(1:3), gammaij(1:6)/)
  !
END SUBROUTINE PDECons2Prim

! Basically a helper...
RECURSIVE SUBROUTINE MatrixInverse3x3(M,iM,det) 
    !---------------
    ! compute the determinant det of the NxN-matrix M with N=3.
    !---------------
    IMPLICIT NONE
    ! input variables 
    REAL, INTENT(IN)   :: M(3,3)
    ! output variables
    REAL, INTENT(OUT)    :: iM(3,3)
    REAL, INTENT(OUT)    :: det
    ! output variables
    REAL    :: Id(3,3)
    INTEGER :: i,j
    ! 
    det = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-M(2,1)*M(1,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*M(1,2)*M(2,3)-M(3,1)*M(1,3)*M(2,2)
    IF(det*det.LT.1e-20) THEN
        print *, 'FATAL ERROR: det = 0'
        CALL EXAHYPE_ABORT
    ENDIF
    !
    iM(1,1) =M(2,2)*M(3,3)-M(2,3)*M(3,2)
    iM(1,2) =M(1,3)*M(3,2)-M(1,2)*M(3,3)
    iM(1,3) =M(1,2)*M(2,3)-M(1,3)*M(2,2)
    iM(2,1) =M(2,3)*M(3,1)-M(2,1)*M(3,3)
    iM(2,2) =M(1,1)*M(3,3)-M(1,3)*M(3,1)
    iM(2,3) =M(1,3)*M(2,1)-M(1,1)*M(2,3)
    iM(3,1) =M(2,1)*M(3,2)-M(2,2)*M(3,1)
    iM(3,2) =M(1,2)*M(3,1)-M(1,1)*M(3,2)
    iM(3,3) =M(1,1)*M(2,2)-M(1,2)*M(2,1)
    iM = iM/det
    !
    Id = MATMUL(M,iM)
    DO i=1,3
        DO j=1,3
            IF(i.eq.j) THEN
                IF((Id(i,j)-1.)**2..GT.1e-18) THEN
                    print *, 'FATAL ERROR 2: det = 0'
                    CALL EXAHYPE_ABORT
                ENDIF
            ELSE
                IF((Id(i,j)**2).GT.1e-18) THEN
                    print *, 'FATAL ERROR 3: det = 0'
                    CALL EXAHYPE_ABORT
                ENDIF
            ENDIF
        ENDDO
    ENDDO
    !
    CONTINUE
    !
END SUBROUTINE MatrixInverse3x3
