/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifndef _EXAHYPE_SOLVERS_TEMPORARY_VARIABLES_H_
#define _EXAHYPE_SOLVERS_TEMPORARY_VARIABLES_H_

#include <vector>

namespace exahype {
  namespace solvers {
    class PredictionTemporaryVariables;
    class MergingTemporaryVariables;

    //helper function to allocate and deallocate memory
    double* allocateArray( std::vector<int>& heapIndices, const int size );
    void freeArrays( std::vector<int>& heapIndices );
    
    /**
     * Initialises temporary variables
     * used for the Prediction and SolutionRecomputation
     * mapping.
     *
     * \note We parallelise over the domain
     * (mapping is copied for each thread) and
     * over the solvers registered on a cell.
     *
     * \note We need to initialise the temporary variables
     * per mapping and not in the solvers since the
     * solvers in exahype::solvers::RegisteredSolvers
     * are not copied for every thread.
     */
    void initialiseTemporaryVariables(PredictionTemporaryVariables& temporaryVariables);

    /**
     * Deletes temporary variables
     * used for the Prediction and SolutionRecomputation
     * mapping.
     */
    void deleteTemporaryVariables(PredictionTemporaryVariables& temporaryVariables);

    /**
     * Initialises temporary variables
     * used in the Merging and SolutionRecomputation
     * mapping.
     *
     * \note We parallelise over the domain
     * (mapping is copied for each thread) and
     * over the solvers registered on a cell.
     *
     * \note We need to initialise the temporary variables
     * per mapping and not in the solvers since the
     * solvers in exahype::solvers::RegisteredSolvers
     * are not copied for every thread.
     */
    void initialiseTemporaryVariables(MergingTemporaryVariables& temporaryVariables);

    /**
     * Deletes temporary variables
     * used in the Merging and SolutionRecomputation
     * mapping.
     */
    void deleteTemporaryVariables(MergingTemporaryVariables& temporaryVariables);
  }
}

struct exahype::solvers::PredictionTemporaryVariables { // TODO(Dominic): Realise per solver.
private:
    /**
     * This is a container, so it is not a good idea to copy it.
     */
  PredictionTemporaryVariables( const PredictionTemporaryVariables& ) {};
public:
  PredictionTemporaryVariables();

  std::vector<int> _dataHeapIndices;

  /**
   * Per solver, temporary variables for storing degrees of freedom of space-time predictor
   * sized variables.
   */
  double*** _tempSpaceTimeUnknowns = nullptr;
  /**
   * Per solver, temporary variables for storing degrees of freedom of space-time
   * volume flux sized variables.
   */
  double*** _tempSpaceTimeFluxUnknowns = nullptr;

  /**
   * Per solver, temporary variables for storing degrees of freedom of solution
   * sized variables.
   *  // TODO(Dominic): This variable can be eliminated from the nonlinear kernels.
   */
  double** _tempUnknowns = nullptr;

  /**
   * Per solver, temporary variables for storing degrees of freedom of volume flux
   * sized variables.
   *  // TODO(Dominic): This variable can be eliminated from the nonlinear kernels.
   */
  double** _tempFluxUnknowns = nullptr;

  //TODO KD describe what it is
  double*** _tempPointForceSources = nullptr;
};

struct exahype::solvers::MergingTemporaryVariables {
private:
    /**
     * This is a container, so it is not a good idea to copy it.
     */
  MergingTemporaryVariables( const MergingTemporaryVariables& ) {};
public:
  MergingTemporaryVariables();

  std::vector<int> _dataHeapIndices;

  /**
   * Temporary variable per solver for storing
   * space-time face unknowns.
   */
  //  double**  _tempSpaceTimeFaceUnknownsArray  = nullptr; todo

  /**
   * Temporary variable per solver for storing
   * face unknowns.
   */
  double***  _tempFaceUnknowns = nullptr;
};

#endif
