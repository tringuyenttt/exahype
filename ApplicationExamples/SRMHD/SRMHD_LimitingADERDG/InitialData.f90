! MHD Initial Data

! InitialData is chosen by the C function pointer "idfunc" which is set
! in MHDSolver_ADERDG.cpp in the init() method.

SUBROUTINE InitialJet(x, Q)
    ! Computes the MHD Jet initial data.
    ! This is also used for boundary condition value computation.
    !
    ! The Jet enters the domain at the zy plane (x=0) in a small circle
    ! with radius rb (=1) around the origin.
    ! A small grid setup is for y,z = min:-5, max:+5
    ! A big grid setup is for y,z = min:-25, max:+25
    ! In any case, the extend in x should be large (5 to 25).
    ! Also use long simulation times (like t=40).
    ! 
    ! This simulation was made for 3D but also works in 2D. The jet
    ! then just enters at the y axis around |y|<1.
    ! 
    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    ! Asked for initial data: With this boolean variable we could
    ! determine whether this is a call for setting ID on the whole
    ! domain or only BC on the yz axis.
    !LOGICAL, INTENT(IN)            :: askedForInitialData
    ! Local variables
    REAL :: rhoa, pa, va, Ms, eta
    REAL :: rhob, pb, vb, betab
    REAL :: V(nVar), BV(3)
    REAL :: rho, rb, xlim

    Ms   = 4.0
    eta  = 1.e-2
    betab= 10.0
    
    rhoa = 1.0
    va   = 0.0

    rhob = rhoa*eta
    vb   = 0.99
    
    pa   = eta*abs(vb)**2 / (gamma*(gamma -1.0)*Ms**2  - gamma*abs(vb)**2)     
    
    BV(1) = sqrt(2.0*Pa/betab)
    BV(2) = 0.0
    BV(3) = 0.0

    rho = SQRT(SUM(x(2:nDim)**2)) ! radius on yz plane, or radius in y
    rb = 1.0 ! Beam radius
    xlim = 0.2 ! Some artificial small area above the yz plane

    ! cf the same query at PDE.f90/MHDJetBC
    IF (x(1).LT.xlim.and.rho.LE.rb) THEN
      ! Jet
      V = (/ rhob, vb, 0.0, 0.0, pa, BV(1:3), 0.0 /)      
    ELSE
      ! vacuum
      V = (/ rhoa, va, 0.0, 0.0, pa, BV(1:3), 0.0 /)
    END IF 
    !
    ! Now convert to conservative variables
    !
    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE InitialJet

SUBROUTINE InitialAlfenWave(x, Q)
    ! Get the AlfenWave for t=0
    ! This subroutine is neccessary for a common argument interface InitialFooBar(x,Q)
    ! for all initial data.
    USE, INTRINSIC :: ISO_C_BINDING
    Use Parameters, ONLY : nDim, nVar
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    CALL AlfenWave(x, Q, 0.0)
END SUBROUTINE InitialAlfenWave

  
SUBROUTINE AlfenWave(x, Q, t)
    ! Computes the AlfenWave at a given time t.
    ! Use it ie. with t=0 for initial data
    ! Use it for any other time ie. for comparison

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)

    rho0 = 1.
    p0   = 1.
    eta  = 1.
    B0   = 1.0 
    !
    hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0
    tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
    tempab = 2.0 * eta * B0**2 / tempaa
    tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
    va2 = b0**2 / ( tempaa * tempac)
    vax = sqrt ( va2)
    !
    BV(1) = B0
    BV(2) = eta * B0 * COS(2*Pi*( x(1) - vax*t))
    BV(3) = eta * B0 * SIN(2*Pi*( x(1) - vax*t))
    !
    VV(1)   = 0.0
    VV(2:3) = - vax * BV(2:3) / B0
    !
    ! Now convert to conservative variables
    !
    V = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)
    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE AlfenWave
  

SUBROUTINE InitialBlast(x, Q)
    ! Blast wave: See 10.1137/0915019 or 10.1016/j.cpc.2014.03.018
    ! Simulation domain: 0.0 .. 1.0 (squared)

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL :: xR(nDim), rho0, p0, p1, rho1, r, r0, r1, taper
    REAL :: V(nVAR)

    ! Here, we set the primitive variables 
    V(1) = 0.0   ! rho
    V(2) = 0.0   ! vx
    V(3) = 0.0   ! vy
    V(4) = 0.0   ! vz
    V(5) = 0.0   ! p
    V(6) = 1/sqrt(2.0)   ! bx
    V(7) = 1/sqrt(2.0)   ! by
    V(8) = 0.0   
    V(9) = 0.0   ! preserving

    xR(1) = x(1) - 0.5
    xR(2) = x(2) - 0.5

    p0   = 10
    p1   = 0.1
    rho0 = 1.0
    rho1 = 1.0

    r = SQRT(xR(1)**2 + xR(2)**2)
    r0 = 0.1
    r1 = 1.0

    !
    IF(r .LE. r0) THEN
        V(1) = rho0
        V(5) = p0
!    ELSEIF ( r.GT.r0 .AND. r.LE.1) THEN
!        taper = (1.-(r-r0)/(r1-r0))
!        V(1) = rho0 + taper*(rho1-rho0)
!        taper = (1.-(r-r0)/(r1-r0))
!        V(5) = p0 + taper*(p1-p0)
    ELSE  
        V(1) = rho1
        V(5) = p1
    ENDIF

    CALL PDEPrim2Cons(Q,V)
END SUBROUTINE

SUBROUTINE InitialOrsagTang(x, Q)
    ! Simulation Domain: 0 .. 2*pi = 6.283185307179586

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 

    REAL :: rho0, p0, vel0, B0
    REAL :: V(nVAR), VV(3), BV(3)

    ! RMHDOrszagTang
    rho0 = 1.0
    p0   = 1.0
    vel0 = 0.75
    B0   = 1.0
    !
    VV(1) = -1./SQRT(2.0)*vel0*SIN(x(2))
    VV(2) =  1./SQRT(2.0)*vel0*SIN(x(1))
    VV(3) = 0. 
    !
    BV(1) = -B0*SIN(x(2))
    BV(2) =  B0*SIN(2.0*x(1))
    BV(3) =  0.
    !
    V = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)

    Call PDEPrim2Cons(Q,V)
END SUBROUTINE InitialOrsagTang

SUBROUTINE InitialRotor(x,Q)
    ! see http://adsabs.harvard.edu/abs/2004MSAIS...4...36D
    ! Domain: 0...1.0  (square domain)
    ! gamma: 5/3

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL :: rho, p, r
    REAL :: v0, rho0, rho1, p0, p1
    REAL :: V(nVAR), VV(3), BV(3)
    
    REAL :: x0(2), r0, r1, f
    
    x0 = (/ 0.5, 0.5 /)
    r0 = 0.1
    r1 = 0.115
    
    rho0 = 10.0
    rho1 = 1.0 

    p0 = 0.5
    p1 = 0.5
 
    v0 = 0.995 ! 1.0 corresponds to the speed of light
    r = sqrt ( (x(1)-x0(1))**2 + (x(2)-x0(2))**2 )   ! cylindrical

    IF  (r < r0) THEN
        rho   =  rho0                   
        VV(1) = -v0 * (x(2)-x0(2))/r0 ! angular velocity v_phi = w * r * e_phi , with e_phi = (-y/r, x/r )^T
        VV(2) =  v0 * (x(1)-x0(1))/r0
        VV(3) =  0.
        p     =  p1
!    ELSE IF (r.GE.r0.AND.r.LT.r1) THEN
!        f = (r1 - r) / (r1-r0) ! linear interpolant;
!
!        rho   = rho1 + (rho0-rho1) * f
!        VV(1) = -v0 * f * (x(2)-x0(2))/r
!        VV(2) =  v0 * f * (x(1)-x0(1))/r
!        VV(3) =  0.
!        p     =  p1
    ELSE
        rho   = rho1
        VV    = 0.
        p     = p0
    ENDIF
    !
    BV(1) = 1.0
    BV(2) = 0.
    BV(3) = 0.
    !
    V = (/ rho, VV(1:3), p, BV(1:3), 0.0 /)
    CALL PDEPrim2Cons(Q, V)
END SUBROUTINE InitialRotor

SUBROUTINE InitialShockTube(x,Q)
    ! see Zanotti, O., Dumbser, M., 2015. A high order special relativistic hydrodynamic and magnetohydrodynamic code with space–time adaptive mesh refinement. Computer Physics Communications 188, 110–127.
    ! Domain: 0...1.0  (square domain)
    ! gamma: 5/3

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    
    REAL :: V(nVAR)
    
    IF  (x(1).LT.0.5) THEN
!       V = (rho (vx,vy,vz) p (bx,by,bz), div_cleaning)
       V = (/ 1.08, 0.4,0.3,0.2, 0.95, 2.0,0.3,0.3, 0.0 /)
    ELSE
       V = (/ 1.0, -0.45,0.2,0.2, 1.0, 2.0,-0.7,0.5, 0.0 /)
    ENDIF
    CALL PDEPrim2Cons(Q, V)
END SUBROUTINE InitialShockTube
