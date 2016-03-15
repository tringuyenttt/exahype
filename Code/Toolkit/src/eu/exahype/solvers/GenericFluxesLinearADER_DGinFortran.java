package eu.exahype.solvers;

public class GenericFluxesLinearADER_DGinFortran implements Solver {
  public static final String Identifier = GenericFluxesLinearADER_DGinC.Identifier;

  private int _dimensions;
  private int _numberOfVariables;
  private int _order;

  public GenericFluxesLinearADER_DGinFortran(int dimensions, int numberOfVariables, int order) {
    _dimensions = dimensions;
    _numberOfVariables = numberOfVariables;
    _order = order;
  }
  public void writeHeader(int dimensions, int numberOfVariables, int order) {
    _dimensions = dimensions;
    _numberOfVariables = numberOfVariables;
    _order = order;
  }
  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo
    Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName);

    writer.write("  private:\n");
    if (_dimensions == 2) {
      writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
    } else {
      writer.write(
          "    static void flux(const double* const Q, double* f, double* g, double* h);\n");
    }
    writer.write(
        "    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write(
        "    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");

    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::spaceTimePredictorLinear<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::volumeIntegralLinear( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::surfaceIntegralLinear( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::riemannSolverLinear<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "   return kernels::aderdg::generic::fortran::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    writer.write(
        "   kernels::aderdg::generic::fortran::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(
        solverName, writer, projectName, _numberOfVariables, _order);

    int digits = String.valueOf(_numberOfVariables).length();

    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfVariables) + "\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    int digits = String.valueOf(_numberOfVariables).length();

    writer.write("SUBROUTINE PDEEigenvalues(Lambda,Q,nv) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar), nv(d)  \n");
    writer.write("  REAL, INTENT(OUT) :: Lambda(nVar)  \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  Lambda(" + String.format("%" + digits + "d", i + 1) + ") = 0.0\n");
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEEigenvalues\n");

    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDEFlux(F,Q) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar) \n");
    writer.write("  REAL, INTENT(OUT) :: F(nVar,d) \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    writer.write("  !\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 1) = 0.0\n");
    }
    writer.write("  !\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 2) = 0.0\n");
    }
    if (_dimensions == 3) {
      writer.write("  !\n");
      for (int i = 0; i < _numberOfVariables; i++) {
        writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 3) = 0.0\n");
      }
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEFlux \n");
    
    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDENCP(BgradQ,Q,gradQ) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d) \n");
    writer.write("  REAL, INTENT(OUT) :: BgradQ(nVar,d) \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    writer.write("  !\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  BgradQ(" + String.format("%" + digits + "d", i + 1) + ", 1) = 0.0\n");
    }
    writer.write("  !\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      writer.write("  BgradQ(" + String.format("%" + digits + "d", i + 1) + ", 2) = 0.0\n");
    }
    if (_dimensions == 3) {
      writer.write("  !\n");
      for (int i = 0; i < _numberOfVariables; i++) {
        writer.write("  BgradQ(" + String.format("%" + digits + "d", i + 1) + ", 3) = 0.0\n");
      }
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDENCP \n");    
    
    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDEMatrixB(Bn,Q,nv) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar), nv(d) \n");
    writer.write("  REAL, INTENT(OUT) :: Bn(nVar,nVar)  \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    writer.write("  !\n");
    for (int i = 0; i < _numberOfVariables; i++) {
      for (int j = 0; j < _numberOfVariables; j++) {
        writer.write("  Bn(" + String.format("%" + digits + "d", i + 1) + ", " + String.format("%" + digits + "d", j + 1) + ") = 0.0\n");
      }
      writer.write("  !\n");
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEMatrixB \n");     }
  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement

    writer.write("  MODULE typesDef \n");
    writer.write("    IMPLICIT NONE  \n");
    writer.write("    PUBLIC  \n");
    writer.write("    ! \n");
    writer.write(
        "    ! ================================== This part of the typesDef can be modified by the user.  ==================================  \n");
    writer.write("    ! \n");
    writer.write(
        "    INTEGER, PARAMETER             :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !!  \n");
    writer.write(
        "    INTEGER, PARAMETER             :: N = " + _order + "                               ! Polynomial degree of our approximation in space and time  \n");
    writer.write(
        "    INTEGER, PARAMETER             :: nDim = "+ _dimensions + "                            ! The number of space dimensions that we actually want to simulate  \n");
    writer.write(
        "    INTEGER, PARAMETER             :: nVar = "+ _numberOfVariables + "                            ! The number of variables of the PDE system  \n");
    writer.write(
        "    INTEGER, PARAMETER             :: nDOF(0:3) = (/ " + (_order+1) + ", " + (_order+1) + ", " + (_order+1) + ", " + (_order+1) + " /)                           ! The number of degrees of freedom in space and time  \n");
    writer.write("     \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: wGPN(N+1)     = (/ 0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: xiGPN(N+1)    = (/ 0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: F0(N+1)       = (/ 1.526788125457e+00, -8.136324494869e-01, 4.007615203117e-01, -1.139171962820e-01 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: FLcoeff(N+1)  = (/ 1.52678812545727, -0.813632449486927, 0.400761520311650, -0.113917196281990 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: FRcoeff(N+1)  = (/ -0.113917196281990, 0.400761520311651, -0.813632449486928, 1.52678812545727 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: dudx(N+1,N+1) = reshape( (/ -6.66400047270456, -1.51511522959847, 0.657396448516548, -1.16125633832453, 9.72030883137039, -0.768828784446417, -2.94134046256143, 4.21756469699036, -4.21756469699036, 2.94134046256143, 0.768828784446416, -9.72030883137039, 1.16125633832453, -0.657396448516549, 1.51511522959847, 6.66400047270456 /), (/N+1,N+1 /) ) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: Kxi(N+1,N+1)  = reshape( (/ -1.15905242621428, 1.69062826161229, -0.733550157264387, 0.201974321866383, -0.494037528020548, -0.250693983347796, 0.959090465730300, -0.214358954361956, 0.214358954361956, -0.959090465730300, 0.250693983347796, 0.494037528020548, -0.201974321866383, 0.733550157264387, -1.69062826161229, 1.15905242621428 /), (/N+1,N+1 /) ) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: iK1(N+1,N+1)  = reshape( (/ 0.546435362419645, 1.01885331677130, 1.02401050669309, 0.974005058264396, -0.144326183293257, 0.584759972857323, 1.00074377855320, 1.02401050669309, 0.101462359828863, -0.170263724267844, 0.584759972857323, 1.01885331677130, -6.687578310368468E-002, 0.101462359828862, -0.144326183293257, 0.546435362419644 /), (/N+1,N+1 /) ) \n");
    writer.write("     \n");
    writer.write("   \n");
    writer.write("    TYPE tFace \n");
    writer.write(
        "      DOUBLE PRECISION, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector  \n");
    writer.write(
        "      DOUBLE PRECISION, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector  \n");
    writer.write(
        "      INTEGER          :: Left, Right                         ! pointer to left and right element  \n");
    writer.write(
        "      DOUBLE PRECISION             :: nv(d)                               ! face normal vector  \n");
    writer.write("    END TYPE       \n");
    writer.write("    TYPE(tFace), POINTER :: Face(:)  \n");
    writer.write("  END MODULE typesDef  \n");
  }
}
