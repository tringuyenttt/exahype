! MHD Boundary conditions for specific problems

! Note that the exact boundary conditions are not specified in Fortran as we
! cannot easily use C function pointers (bcfunc as in BoundaryConditions.h)

! The BC mostly go along with specific InitialData.

SUBROUTINE BoundaryOutflow(x,t,dt,faceIndex,nv,fluxIn,stateIn,fluxOut,stateOut)
  USE Parameters, ONLY : nVar, nDim
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list
  REAL, INTENT(IN) :: x(nDim), t, dt, nv(3)
  INTEGER, INTENT(IN) :: faceIndex
  REAL, INTENT(IN) :: fluxIn(nVar,nDim), stateIn(nVar)
  REAL, INTENT(OUT) :: fluxOut(nVar,nDim), stateOut(nVar)


  ! Outgoing BC
  fluxOut = fluxIn
  stateOut = stateIn
END SUBROUTINE BoundaryOutflow

! AlfenWave BC defined in C as they need kernel:: data structures

SUBROUTINE BoundaryJet(x,t,dt,faceIndex,nv,fluxIn,stateIn,fluxOut,stateOut)
  ! Boundary conditions for the MHDJet
  USE Parameters, ONLY : nVar, nDim
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! Argument list
  REAL, INTENT(IN) :: x(nDim), t, dt, nv(nDim)
  INTEGER, INTENT(IN) :: faceIndex
  REAL, INTENT(IN) :: fluxIn(nVar), stateIn(nVar)
  REAL, INTENT(OUT) :: fluxOut(nVar), stateOut(nVar)
  ! Local variables
  REAL :: rho, rb, allFluxOut(nVar,nDim)

  rho = SQRT(SUM(x(2:nDim)**2)) ! radius on yz plane, or radius in y
  rb = 1.0 ! Beam radius
  
  ! faces: 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
  ! ie.    0 x=0   1 x=max, 2 y=0    3 y=max 4 z=0     5 z=max
  IF (faceIndex.eq.0.and.rho<rb) THEN
    CALL InitialJet(x, stateOut, .FALSE.) ! do jet
    CALL PDEFlux(allFluxOut, stateOut)
    fluxOut = MATMUL(allFluxOut, nv) ! Project flux in boundary direction
  ELSE
    ! Outgoing BC
    fluxOut = fluxIn
    stateOut = stateIn
  ENDIF
END SUBROUTINE BoundaryJet

