! SRHD application: Shocktube parameters

  MODULE Parameters 
    use iso_c_binding
    IMPLICIT NONE  
    PUBLIC  

    ! This typesDef.f90 is a relict still from the old Fortran interface in Exahype.
    ! However, the following two variables are needed in the Fortran code. They
    ! should provided by the glue code generator in an appropriate way.

    ! If you modify the SRHD.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !

    INTEGER, PARAMETER             :: nDim = 2                            ! The number of space dimensions that we actually want to simulate  
    INTEGER, PARAMETER             :: nVar = 5                            ! The number of variables of the PDE system  
    
    ! This was moved from TYPE(tEquation) :: EQN from the MainVariables.f90 module by Olindo
    REAL, PARAMETER                :: gamma = 5.0/3.0
    
    ! Add further parameters as from 
    !REAL(C_FLOAT), BIND(C)  :: foobar
  END MODULE Parameters  
