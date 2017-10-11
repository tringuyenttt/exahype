  MODULE typesDef 
    IMPLICIT NONE  
    PUBLIC  
    ! 
    ! ================================== This part of the typesDef can be modified by the user.  ==================================  
    ! 
    INTEGER, PARAMETER             :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !!  
    INTEGER, PARAMETER             :: N = 3                               ! Polynomial degree of our approximation in space and time  
    INTEGER, PARAMETER             :: nDim = 3                            ! The number of space dimensions that we actually want to simulate  
    INTEGER, PARAMETER             :: nVar = 4                            ! The number of variables of the PDE system  
    INTEGER, PARAMETER             :: nDOF(0:3) = (/ 4, 4, 4, 4 /)                           ! The number of degrees of freedom in space and time  
     
    DOUBLE PRECISION, PARAMETER    :: wGPN(N+1)     = (/ 0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273 /) 
    DOUBLE PRECISION, PARAMETER    :: xiGPN(N+1)    = (/ 0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262 /) 
    DOUBLE PRECISION, PARAMETER    :: F0(N+1)       = (/ 1.526788125457e+00, -8.136324494869e-01, 4.007615203117e-01, -1.139171962820e-01 /) 
    DOUBLE PRECISION, PARAMETER    :: FLcoeff(N+1)  = (/ 1.52678812545727, -0.813632449486927, 0.400761520311650, -0.113917196281990 /) 
    DOUBLE PRECISION, PARAMETER    :: FRcoeff(N+1)  = (/ -0.113917196281990, 0.400761520311651, -0.813632449486928, 1.52678812545727 /) 
    DOUBLE PRECISION, PARAMETER    :: dudx(N+1,N+1) = reshape( (/ -6.66400047270456, -1.51511522959847, 0.657396448516548, -1.16125633832453, 9.72030883137039, -0.768828784446417, -2.94134046256143, 4.21756469699036, -4.21756469699036, 2.94134046256143, 0.768828784446416, -9.72030883137039, 1.16125633832453, -0.657396448516549, 1.51511522959847, 6.66400047270456 /), (/N+1,N+1 /) ) 
    DOUBLE PRECISION, PARAMETER    :: Kxi(N+1,N+1)  = reshape( (/ -1.15905242621428, 1.69062826161229, -0.733550157264387, 0.201974321866383, -0.494037528020548, -0.250693983347796, 0.959090465730300, -0.214358954361956, 0.214358954361956, -0.959090465730300, 0.250693983347796, 0.494037528020548, -0.201974321866383, 0.733550157264387, -1.69062826161229, 1.15905242621428 /), (/N+1,N+1 /) ) 
    DOUBLE PRECISION, PARAMETER    :: iK1(N+1,N+1)  = reshape( (/ 0.546435362419645, 1.01885331677130, 1.02401050669309, 0.974005058264396, -0.144326183293257, 0.584759972857323, 1.00074377855320, 1.02401050669309, 0.101462359828863, -0.170263724267844, 0.584759972857323, 1.01885331677130, -6.687578310368468E-002, 0.101462359828862, -0.144326183293257, 0.546435362419644 /), (/N+1,N+1 /) ) 
     
   
    TYPE tFace 
      DOUBLE PRECISION, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector  
      DOUBLE PRECISION, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector  
      INTEGER          :: Left, Right                         ! pointer to left and right element  
      DOUBLE PRECISION             :: nv(d)                               ! face normal vector  
    END TYPE       
    TYPE(tFace), POINTER :: Face(:)  
  END MODULE typesDef  
