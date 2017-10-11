SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
  USE typesDef, ONLY : nVar, d
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), nv(d) 
  REAL, INTENT(OUT) :: Lambda(nVar) 
  ! Local variables 
  REAL :: c = 1.0
  !
  Lambda = (/ c /)                        ! The eigenvalues of the Euler equations 
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
  REAL :: p, irho, lam, mu 
  REAL :: Qx(nVar), Qy(nVar), Qz(nVar)
  
  BgradQ(:,:) = 0.0
  
  !BgradQ(:,1) = -gradQ(:,1) 
  !BgradQ(:,2) = -gradQ(:,2)
  !BgradQ(:,3) = -gradQ(:,3)

  Qx(:) = gradQ(:,1) 
  Qy(:) = gradQ(:,2)
  Qz(:) = gradQ(:,3)

  BgradQ(1,1) = -Qx(1)
  !BgradQ(2,1) = -Qx(1)

  BgradQ(1,2) = -Qy(1)
  !BgradQ(3,2) = -Qy(1)

  BgradQ(1,3) = -Qz(1)
  !BgradQ(4,3) = -Qz(1) 
  
  !PRINT *, 'PDENCP is not used for nonlinear-solvers' 
  !CALL EXIT  
  !
END SUBROUTINE PDENCP 
 


SUBROUTINE PDEMatrixB(Bn,Q,nv)
  USE typesDef, ONLY : nVar, d
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), nv(d)   
  REAL, INTENT(OUT) :: Bn(nVar,nVar) 
  ! Local variables 
  REAL :: p, irho, lam, mu 
  REAL :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)

  B1 = 0. 
  B2 = 0. 
  B3 = 0.

  B1(1,1) = -1.0
  !B1(2,1) = -1.0

  B2(1,1) = -1.0
  !B2(3,2) = -1.0

  B3(1,1) = -1.0
  !B3(4,3) = -1.0

  Bn = B1*nv(1) + B2*nv(2) + B3*nv(3)
  !
 
  !
END SUBROUTINE PDEMatrixB 
