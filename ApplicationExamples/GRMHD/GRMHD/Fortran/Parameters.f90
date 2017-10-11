! MHD Parameters: Riemann problem, Rotor problem, Blast wave
  
  MODULE Parameters 
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

    ! Here, we obtain DIMENSIONS from the C++ preprocessor
#if defined(Dim3)
    INTEGER, PARAMETER             :: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             :: nDim = 2                   ! The number of space dimensions
#endif
    INTEGER, PARAMETER             :: nVar = 19                           ! The number of variables of the PDE system  
    
    
    ! Ideal EOS:
    ! 4/3 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
    ! 2.0 used for TOV stars
    REAL, PARAMETER                :: gamma = 2.0 ! 4.0/3.0
      
    ! Divergence cleaning:
    ! 1.0 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
    ! Note that currently we don't add the contribution to the source.
    REAL, PARAMETER :: DivCleaning_a = 1.0

  END MODULE Parameters  
