package eu.exahype.solvers;

import eu.exahype.IOUtils;

public class GenericADERDG implements Solver {
  public static final String Identifier = "generic::fluxes";

  private int _dimensions;
  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _order;
//  private int _patchSize;
  private boolean _enableProfiler;
  private boolean _hasConstants;
  private boolean _isLinear;
  private boolean _isFortran;

  public GenericADERDG(int dimensions, int numberOfVariables, int numberOfParameters,
      int order, boolean enableProfiler, boolean hasConstants, boolean isLinear, boolean isFortran) {
    _dimensions         = dimensions;
    _numberOfVariables  = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _order              = order;
//    _patchSize = patchSize;
    _enableProfiler     = enableProfiler;
    _hasConstants       = hasConstants;
    _isLinear           = isLinear;
    _isFortran          = isFortran;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericADERDGSolverHeader.template");
  
	  content = content.replaceAll("\\$\\{Project\\}", projectName);
	  content = content.replaceAll("\\$\\{Solver\\}", solverName);
	  
	  String profilerInclude                     = "";
	  String solverConstructorSignatureExtension = "";
	  if (_enableProfiler) {
		  profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
		  solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler"; 
	  }
	  if (_hasConstants) {
		  solverConstructorSignatureExtension += ", exahype::Parser::ParserView constants"; // TODO(Dominic): Why pass by value? 
	  }
	  content = content.replaceAll("\\$\\{ProfilerInclude\\}",profilerInclude);
	  content = content.replaceAll("\\$\\{SolverConstructorSignatureExtension\\}", solverConstructorSignatureExtension);
	  
	  writer.write(content);
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericADERDGSolverInCGeneratedCode.template");
    
	  content = content.replaceAll("\\$\\{Project\\}", projectName);
	  content = content.replaceAll("\\$\\{Solver\\}", solverName);
	  //
	  String profilerInclude                     = "";
	  String solverConstructorSignatureExtension = "";
	  String solverConstructorArgumentExtension  = "";
	  if (_enableProfiler) {
		  profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
		  solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler";
		  solverConstructorArgumentExtension  += ", std::move(profiler)";
		  
      content = content.replaceAll("\\$\\{BeforeSpaceTimePredictor\\}", "  _profiler->start(\"spaceTimePredictor\");");  
      content = content.replaceAll("\\$\\{AfterSpaceTimePredictor\\}", "  _profiler->stop(\"spaceTimePredictor\");"); 
      content = content.replaceAll("\\$\\{BeforeSolutionUpdate\\}", "  _profiler->start(\"solutionUpdate\");"); 
      content = content.replaceAll("\\$\\{AfterSolutionUpdate\\}", "  _profiler->stop(\"solutionUpdate\");"); 
      content = content.replaceAll("\\$\\{BeforeVolumeIntegral\\}", "  _profiler->start(\"volumeIntegral\");"); 
      content = content.replaceAll("\\$\\{AfterVolumeIntegral\\}", "  _profiler->stop(\"volumeIntegral\");"); 
      content = content.replaceAll("\\$\\{BeforeSurfaceIntegral\\}", "  _profiler->start(\"surfaceIntegral\");"); 
      content = content.replaceAll("\\$\\{AfterSurfaceIntegral\\}", "  _profiler->stop(\"surfaceIntegral\");"); 
      content = content.replaceAll("\\$\\{BeforeRiemannSolver\\}", "  _profiler->start(\"riemannSolver\");"); 
      content = content.replaceAll("\\$\\{AfterRiemannSolver\\}", "  _profiler->stop(\"riemannSolver\");"); 
      content = content.replaceAll("\\$\\{BeforeBoundaryConditions\\}", "  _profiler->start(\"boundaryConditions\");"); 
      content = content.replaceAll("\\$\\{AfterBoundaryConditions\\}", "  _profiler->stop(\"boundaryConditions\");"); 
      content = content.replaceAll("\\$\\{BeforeStableTimeStepSize\\}", "  _profiler->start(\"stableTimeStepSize\");"); 
      content = content.replaceAll("\\$\\{AfterStableTimeStepSize\\}", "  _profiler->stop(\"stableTimeStepSize\");"); 
      content = content.replaceAll("\\$\\{BeforeSolutionAdjustment\\}", "  _profiler->start(\"solutionAdjustment\");"); 
      content = content.replaceAll("\\$\\{AfterSolutionAdjustment\\}", "  _profiler->stop(\"solutionAdjustment\");"); 
      content = content.replaceAll("\\$\\{BeforeFaceUnknownsProlongation\\}", "  _profiler->start(\"faceUnknownsProlongation\");"); 
      content = content.replaceAll("\\$\\{AfterFaceUnknownsProlongation\\}", "  _profiler->stop(\"faceUnknownsProlongation\");"); 
      content = content.replaceAll("\\$\\{BeforeFaceUnknownsRestriction\\}", "  _profiler->start(\"faceUnknownsRestriction\");"); 
      content = content.replaceAll("\\$\\{AfterFaceUnknownsRestriction\\}", "  _profiler->stop(\"faceUnknownsRestriction\");"); 
      content = content.replaceAll("\\$\\{BeforeVolumeUnknownsProlongation\\}", "  _profiler->start(\"volumeUnknownsProlongation\");"); 
      content = content.replaceAll("\\$\\{AfterVolumeUnknownsProlongation\\}", "  _profiler->stop(\"volumeUnknownsProlongation\");"); 
      content = content.replaceAll("\\$\\{BeforeVolumeUnknownsRestriction\\}", "  _profiler->start(\"volumeUnknownsRestriction\");"); 
      content = content.replaceAll("\\$\\{AfterVolumeUnknownsRestriction\\}", "  _profiler->stop(\"volumeUnknownsRestriction\");");
	  } else {
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeSpaceTimePredictor\\}", "");  
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterSpaceTimePredictor\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeSolutionUpdate\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterSolutionUpdate\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeVolumeIntegral\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterVolumeIntegral\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeSurfaceIntegral\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterSurfaceIntegral\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeRiemannSolver\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterRiemannSolver\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeBoundaryConditions\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterBoundaryConditions\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeStableTimeStepSize\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterStableTimeStepSize\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeSolutionAdjustment\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterSolutionAdjustment\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeFaceUnknownsProlongation\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterFaceUnknownsProlongation\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeFaceUnknownsRestriction\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterFaceUnknownsRestriction\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeVolumeUnknownsProlongation\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterVolumeUnknownsProlongation\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{BeforeVolumeUnknownsRestriction\\}", ""); 
      content = content.replaceAll("(\\n|\\r)+\\$\\{AfterVolumeUnknownsRestriction\\}", "");
	  }
	  if (_hasConstants) {
		  solverConstructorSignatureExtension += ", exahype::Parser::ParserView constants"; // TODO(Dominic): Why pass by value? 
	  }
	  
	  content = content.replaceAll("\\$\\{ProfilerInclude\\}",profilerInclude);
	  content = content.replaceAll("\\$\\{SolverConstructorSignatureExtension\\}", solverConstructorSignatureExtension);
	  content = content.replaceAll("\\$\\{SolverConstructorArgumentExtension\\}", solverConstructorArgumentExtension);
	  
	  if (_isLinear) {
	    content = content.replaceAll("\\$\\{NonlinearOrLinear\\}","Linear");
	  } else {
	    content = content.replaceAll("\\$\\{NonlinearOrLinear\\}","Nonlinear");
	  }
	  
	  if (_isFortran) {
      content = content.replaceAll("\\$\\{Language\\}","fortran");
    } else {
      content = content.replaceAll("\\$\\{Language\\}","c");
    }
	  
	  content = content.replaceAll("\\$\\{NumberOfVariables\\}", String.valueOf(_numberOfVariables));
	  content = content.replaceAll("\\$\\{NumberOfParameters\\}",String.valueOf( _numberOfParameters));
	  content = content.replaceAll("\\$\\{Order\\}", String.valueOf(_order));
	  
	  writer.write(content);
  }
  
  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericADERDGSolverInCUserCode.template");
    
    content = content.replaceAll("\\$\\{Project\\}", projectName);
    content = content.replaceAll("\\$\\{Solver\\}", solverName);
    
    content = content.replaceAll("\\$\\{Elements\\}",  String.valueOf( _numberOfParameters+_numberOfVariables));
    content = content.replaceAll("\\$\\{Dimensions\\}",String.valueOf(_dimensions));
    
    int digits = String.valueOf(_numberOfVariables + _numberOfParameters).length();

    String adjustedSolutionValues = "  // State variables:\n";
    for (int i = 0; i < _numberOfVariables; i++) {
      adjustedSolutionValues += "  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) adjustedSolutionValues += "\n";
    }
    if (_numberOfParameters>0) {
      adjustedSolutionValues += "  // Material parameters:\n";
      for (int i = 0; i < _numberOfParameters; i++) {
        adjustedSolutionValues += "  Q[" + String.format("%" + digits + "d", _numberOfVariables+i) + "] = 0.0;";
        if (i<_numberOfParameters-1) adjustedSolutionValues += "\n";
      }
    }

    String eigenvalues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      eigenvalues += "  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) eigenvalues += "\n";
    }

    String flux = "";
    for (int d=0; d<_dimensions; ++d) {
      for (int i = 0; i < _numberOfVariables; i++) {
        flux += "  F["+d+"][" + String.format("%" + digits + "d", i) + "] = 0.0;";
        if (i<_numberOfVariables-1) flux += "\n";
      }
      if (d<_dimensions-1) {
        flux += "\n\n";    
      }
    }
    
    String source = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      source += "  S[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) source += "\n";
    }
    
    String boundaryValues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      boundaryValues += "  stateOut[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) boundaryValues += "\n";
    }
    boundaryValues += "\n\n";
    for (int i = 0; i < _numberOfVariables; i++) {
      boundaryValues += "  fluxOut[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) boundaryValues += "\n";
    }
    
    String ncp = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      ncp += "  BgradQ[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) ncp += "\n";
    }
    
    String matrixb = "";
    for (int i = 0; i < _numberOfVariables*_numberOfVariables; i++) {
      matrixb += "  Bn[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables*_numberOfVariables-1) matrixb += "\n";
    }
    
    content = content.replaceAll("\\$\\{AdjustedSolutionValues\\}",adjustedSolutionValues);
    content = content.replaceAll("\\$\\{Eigenvalues\\}",eigenvalues);
    content = content.replaceAll("\\$\\{Flux\\}",flux);
    content = content.replaceAll("\\$\\{Source\\}",source);
    content = content.replaceAll("\\$\\{BoundaryValues\\}",boundaryValues);
    content = content.replaceAll("\\$\\{NonConservativeProduct\\}",ncp);
    content = content.replaceAll("\\$\\{MatrixB\\}",matrixb);
    
    writer.write(content);
  }

  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a PDF.f90.\n");
  }


  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }
}