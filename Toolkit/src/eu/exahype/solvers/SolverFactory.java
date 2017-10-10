package eu.exahype.solvers;

import java.util.Set;
import java.util.List;
import java.util.Arrays;

import eu.exahype.kernel.ADERDGKernel;
import eu.exahype.kernel.FiniteVolumesKernel;

/**
 * Creates a solver
 * 
 * The idea is that you pass in a kernel object (either FiniteVolumesKernel or 
 * ADERDGKernel) which returns an enum describing the solver type. This 
 * information together with additional properties then is used to instantiate
 * the correct solver class.
 * 
 * <h2> Add support for a new solver type </h2>
 *
 * Switch to either the FiniteVolumesKernel or ADERDGKernel. Add your new solver 
 * type to the enum called KernelType. Next, add a new branch to the routine
 * getKernelType() which returns your brand new solver. Finally, extend the 
 * switch statements in this class to instantiate this new solver.
 *
 */
public class SolverFactory {
  private String _projectName;
  private int _dimensions;
  private boolean _enableProfiler;
  private boolean _enableDeepProfiler;
  private String _microarchitecture;
  private String _pathToLibxsmm;

  
  public SolverFactory(
      String projectName,
      int dimensions,
      boolean enableProfiler,
      boolean enableDeepProfiler,
      String microarchitecture) {
    _projectName = projectName;  
    _dimensions = dimensions;
    _enableProfiler = enableProfiler;
    _enableDeepProfiler = enableDeepProfiler;
    _microarchitecture = microarchitecture;
  }
  
  /**
   * Generates the writer for an ADER-DG solver
   * 
   * Consult class documentation for further details.
   */
  public Solver createADERDGSolver(String solvername, ADERDGKernel kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int order,boolean hasConstants) {
    try { //some solver initialisation can throw IllegalArgumentException if the options are wrong or IOException
      switch (kernel.getKernelType()) {
        case GenericADERDG: 
          return new eu.exahype.solvers.GenericADERDG(_projectName, solvername, _dimensions,
            numberOfVariables, numberOfParameters, namingSchemeNames, order, _enableProfiler, hasConstants, isFortran, kernel );
        case OptimisedADERDG:
          return new eu.exahype.solvers.OptimisedADERDG(_projectName, solvername, _dimensions,
            numberOfVariables, numberOfParameters, namingSchemeNames, order, _microarchitecture,
            _enableProfiler, _enableDeepProfiler, hasConstants, kernel);
      }
    	
/*      
      // TODO JMG Clean
      // else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )) {
        // return new eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC(_dimensions,
            // numberOfVariables, numberOfParameters, order, _microarchitecture,
            // _enableProfiler, hasConstants);
      // }
      else if (!isFortran && kernel.isKernelType( eu.exahype.solvers.OptimisedADERDG.Identifier )) {
        return new eu.exahype.solvers.OptimisedADERDG(_projectName, solvername, _dimensions,
            numberOfVariables, numberOfParameters, namingSchemeNames, order, _microarchitecture,
            _enableProfiler, _enableDeepProfiler, hasConstants, kernel);
        
      }
      // TODO JMG Clean when confirmed legacy
      // else if (!isFortran && kernel.equals( eu.exahype.solvers.KernelEuler2d.Identifier )) {
        // return new eu.exahype.solvers.KernelEuler2d(_projectName, solvername);
      // }
      else if (isFortran && kernel.isKernelType( eu.exahype.solvers.UserDefinedADER_DGinFortran.Identifier )) {
        return new eu.exahype.solvers.UserDefinedADER_DGinFortran(_projectName, solvername);
      }
      else if (!isFortran && kernel.isKernelType( eu.exahype.solvers.UserDefinedADER_DGinC.Identifier )) {
        return new eu.exahype.solvers.UserDefinedADER_DGinC(_projectName, solvername, numberOfVariables,
            numberOfParameters, order, hasConstants, _enableProfiler);
      }
*/      
      System.err.println("ERROR: solver configuration is not supported: "+kernel.toString() );
      return null;
    } catch(Exception e) {
      System.err.println("ERROR: can't create the solver. Error: "+e );
      return null;
    }
  }
  
  /**
   * Generates the writer for an Finite Volumes solver
   * 
   * Consult class documentation for further details.
   */
  public Solver createFiniteVolumesSolver(String solvername, FiniteVolumesKernel kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int patchSize,boolean hasConstants) {
    try { //some solver initialisation can throw IllegalArgumentException if the options are wrong or IOException
      switch (kernel.getKernelType()) {
        case GenericMUSCLHancock: 
          if (isFortran) {
        	// @todo Does not exist yet
            //return new eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInFortran(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
          }
          else {
            return new eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants, kernel);
          }
          break;
        case GenericGodunov: 
            if (isFortran) {
              // @todo Does not exist yet
              //return new eu.exahype.solvers.GenericFiniteVolumesGodunovInFortran(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
            }
            else {
              return new eu.exahype.solvers.GenericFiniteVolumesGodunovInC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants, kernel);
            }
            break;
        case UserDefined: 
            if (isFortran) {
              return new eu.exahype.solvers.UserDefinedFiniteVolumesinFortran(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
            }
            else {
              return new eu.exahype.solvers.UserDefinedFiniteVolumesinC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
            }
        }
      System.err.println("ERROR: solver configuration is not supported: "+kernel.toString() );
      return null;
    } catch(Exception e) {
      System.err.println("ERROR: can't create the solver. Error: "+e );
      return null;
    }
  }
}
