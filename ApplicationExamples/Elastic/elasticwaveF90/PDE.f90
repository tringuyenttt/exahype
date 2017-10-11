SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
  USE typesDef, ONLY : nVar, d
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), nv(d) 
  REAL, INTENT(OUT) :: Lambda(nVar) 
  ! Local variables
  REAL :: lam, mu, rho
  REAL :: cs, cp

  rho = 1.0
  mu = 1.0
  lam = 2.0
  !
  cs = sqrt(mu/rho)
  cp = sqrt((lam+2.0*mu)/rho)
  
  Lambda = (/-cp, -cs, -cs, 0.0, 0.0, 0.0, cp, cs, cs/)                        ! The eigenvalues of the Euler equations 
  !
END SUBROUTINE PDEEigenvalues


SUBROUTINE PDEFlux(F,Q)
   ! not used  
  PRINT *, 'PDEFlux is not used for linear-solvers' 
  CALL EXIT  
 
  !
END SUBROUTINE PDEFlux 
 


SUBROUTINE PDENCP(BgradQ,Q,gradQ)
  USE typesDef, ONLY : nVar, d 
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d)  
  REAL, INTENT(OUT) :: BgradQ(nVar,d) 
  ! Local variables 
  REAL :: rho, lam, mu 
  REAL :: Qx(nVar), Qy(nVar), Qz(nVar)

  rho = 1.0
  mu = 1.0
  lam = 2.0

  !(vx, vy, vz, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
  !  1  2   3    4     5     6     7     8     9
  
  BgradQ(:,:) = 0.0
  
  !BgradQ(:,1) = -gradQ(:,1) 
  !BgradQ(:,2) = -gradQ(:,2)
  !BgradQ(:,3) = -gradQ(:,3)

  Qx(:) = gradQ(:,1) 
  Qy(:) = gradQ(:,2)
  Qz(:) = gradQ(:,3)

  BgradQ(1,1) = -1.0/rho*Qx(4)
  BgradQ(2,1) = -1.0/rho*Qx(7)
  BgradQ(3,1) = -1.0/rho*Qx(8)


  BgradQ(4,1) = -(2.0*mu + lam)*Qx(1)
  BgradQ(5,1) = -lam*Qx(1)
  BgradQ(6,1) = -lam*Qx(1)

  BgradQ(7,1) = -mu*Qx(2)
  BgradQ(8,1) = -mu*Qx(3)



  BgradQ(1,2) = -1.0/rho*Qy(7)
  BgradQ(2,2) = -1.0/rho*Qy(5)
  BgradQ(3,2) = -1.0/rho*Qy(9)


  BgradQ(4,2) = -lam*Qy(2) 
  BgradQ(5,2) = -(2.0*mu + lam)*Qy(2)
  BgradQ(6,2) = -lam*Qy(2)

  BgradQ(7,2) = -mu*Qy(1)
  BgradQ(9,2) = -mu*Qy(3)


  BgradQ(1,3) = -1.0/rho*Qz(8)
  BgradQ(2,3) = -1.0/rho*Qz(9)
  BgradQ(3,3) = -1.0/rho*Qz(6)


  BgradQ(4,3) = -lam*Qz(3) 
  BgradQ(5,3) = -lam*Qz(3)
  BgradQ(6,3) = -(2.0*mu + lam)*Qz(3)

  BgradQ(8,3) = -mu*Qz(1)
  BgradQ(9,3) = -mu*Qz(2)
  

  !CALL EXIT  
  
END SUBROUTINE PDENCP 
 


SUBROUTINE PDEMatrixB(Bn,Q,nv)
  USE typesDef, ONLY : nVar, d
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), nv(d)   
  REAL, INTENT(OUT) :: Bn(nVar,nVar) 
  ! Local variables 
  REAL :: rho, lam, mu 
  REAL :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)
   
  rho = 1.0
  mu =  1.0
  lam = 2.0

  B1 = 0. 
  B2 = 0. 
  B3 = 0.

  B1(1,4) = -1.0/rho
  B1(2,7) = -1.0/rho
  B1(3,8) = -1.0/rho

  B1(4,1) = -(2.0*mu + lam)
  B1(5,1) = -lam
  B1(6,1) = -lam

  B1(7,2) = -mu
  B1(8,3) = -mu


  B2(1,7) = -1.0/rho
  B2(2,5) = -1.0/rho
  B2(3,9) = -1.0/rho

  B2(4,2) = -lam
  B2(5,2) = -(2.0*mu + lam)
  B2(6,2) = -lam

  B2(7,1) = -mu
  B2(9,3) = -mu


  B3(1,8) = -1.0/rho
  B3(2,9) = -1.0/rho
  B3(3,6) = -1.0/rho


  B3(4,3) = -lam 
  B3(5,3) = -lam
  B3(6,3) = -(2.0*mu + lam)

  B3(8,1) = -mu
  B3(9,2) = -mu

 

  Bn = B1*nv(1) + B2*nv(2) + B3*nv(3)

  !print *, Bn, nv(1),  nv(2),  nv(3)
  !
  
  !
END SUBROUTINE PDEMatrixB 
