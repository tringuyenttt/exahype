SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
  USE typesDef, ONLY : nVar, d
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar), nv(d) 
  REAL, INTENT(OUT) :: Lambda(nVar) 
  ! Local variables 
  REAL :: p, u, c 
  !
  u = ( Q(2)*nv(1) + Q(3)*nv(2) + Q(4)*nv(3) )/Q(1)       ! normal velocity 
  !p = (EQNgamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
  p  = (0.4       )*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
  !c = SQRT(EQNgamma*p/Q(1))                              ! sound speed
  c  = SQRT(1.4     *p/Q(1))                              ! sound speed
  !
  Lambda = (/ u-c, u, u, u, u+c /)                        ! The eigenvalues of the Euler equations 
  !
END SUBROUTINE PDEEigenvalues


SUBROUTINE PDEFlux(F,Q)
  USE typesDef, ONLY : nVar, d
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list 
  REAL, INTENT(IN)  :: Q(nVar) 
  REAL, INTENT(OUT) :: F(nVar,d) 
  ! Local variables 
  REAL :: p, irho  
  !
  ! 3D compressible Euler equations 
  !
  irho = 1.0/Q(1)
  !p = (EQNgamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
  p  = (0.4       )*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
  ! 
  F(1,1) = Q(2) 
  F(2,1) = irho*Q(2)*Q(2) + p 
  F(3,1) = irho*Q(2)*Q(3)
  F(4,1) = irho*Q(2)*Q(4)
  F(5,1) = irho*Q(2)*(Q(5)+p)  
  !
  F(1,2) = Q(3) 
  F(2,2) = irho*Q(3)*Q(2)  
  F(3,2) = irho*Q(3)*Q(3) + p 
  F(4,2) = irho*Q(3)*Q(4)
  F(5,2) = irho*Q(3)*(Q(5)+p)  
  ! 
  F(1,3) = Q(4) 
  F(2,3) = irho*Q(4)*Q(2)  
  F(3,3) = irho*Q(4)*Q(3)  
  F(4,3) = irho*Q(4)*Q(4) + p
  F(5,3) = irho*Q(4)*(Q(5)+p)  
  !
END SUBROUTINE PDEFlux 
 


SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
  ! not used  
  PRINT *, 'PDENCP is not used for nonlinear-solvers' 
  CALL EXIT  
  !
END SUBROUTINE PDENCP 
 


SUBROUTINE PDEMatrixB(Bn,Q,nv) 
  ! not used  
  PRINT *, 'PDEMatrixB is not used for nonlinear-solvers' 
  CALL EXIT  
  !
END SUBROUTINE PDEMatrixB 
