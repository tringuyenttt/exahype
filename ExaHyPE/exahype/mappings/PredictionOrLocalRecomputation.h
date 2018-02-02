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
 
#ifndef EXAHYPE_MAPPINGS_PredictionOrLocalRecomputation_H_
#define EXAHYPE_MAPPINGS_PredictionOrLocalRecomputation_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

#include "exahype/solvers/Solver.h"

#include "exahype/solvers/TemporaryVariables.h"

namespace exahype {
namespace mappings {
class PredictionOrLocalRecomputation;
}
}

/**
 * This mapping is one of three mappings (+adapters) which together perform the limiter
 * status spreading and the local recomputation of troubled cells (and their direct) neighbours.
 *
 * |Mapping                 | Event                  | Action                                            |
 * -------------------------------------------------------------------------------------------------------
 * | LimiterStatusSpreading | touchVertexFirstTime   | Merge the face-wise limiter status between local neighbours.|
 * |                        | enterCell              | Determine a unified value of the merged face-wise limiter status values and write it to every face. (Do not update the cell-wise limiter status.)|
 * |                        | prepareSendToNeighbour | Send the unified ace-wise limiter status to neighbouring ranks.|
 * -------------------------------------------------------------------------------------------------------
 * | LocalRollback          | mergeWithNeighbour     | Receive the neighbour ranks' cells' merged limiter status. Directly update the unified merged limiter status value of cells at the remote boundary |
 * |                        | touchVertexFirstTime   | Merge the limiter status between local neighbours (again).
 * |                        | enterCell              | Determine a unified value of the merged face-wise limiter status values and write it to every face. |
 * |                        |                        | (Do not update the cell-wise limiter status.)                                                       |
 * |                        |                        | Rollback the solution in troubled cells and their next two neighbours.
 * |                        | prepareSendToNeighbour | Based on the unified face-wise limiter status send interface values to the neighbours. |
 * ------------------------------------------------------------------------------------------------------
 * |PredictionOrLocalRecomputation | mergeWithNeighbour     | Based on the unified face-wise limiter status receive or drop the interface values send by the neighbours.
 * |                               | touchVertexFirstTIme   | Based on the unified face-wise limiter status merge local direct neighbours.
 * |                               | enterCell              | Recompute the solution in the troubled cells (and their direct neighbours).
 * |                               |                        | Set the cell-wise limiter status to the unified face-wise limiter status.
 * |                               |                        | (The normal time marching does only consider the cell-wise limiter status from now on
 * |                               |                        | and overwrites the face-wise values uniformly with Ok or Troubled after
 * |                               |                        | evaluating the discrete maximum principle (DMP) and the physical admissibility detection (PAD).)
 * ------------------------------------------------------------------------
 *
 * @author Dominic Charrier
 */
class exahype::mappings::PredictionOrLocalRecomputation {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  #ifdef Debug // TODO(Dominic): Exclude shared memory etc.
  /*
   *  Counter for the interior face solves for debugging purposes.
   */
  int _interiorFaceMerges;
  /*
   *  Counter for the boundary face solves for debugging purposes.
   */
  int _boundaryFaceMerges;
  #endif

  /**
   * Local copy of the state which
   * is used to determine if a solver
   * is active in the current algorithm section.
   * (See exahype::runners::Runner for locations
   * where the algorithm section is set. The new
   * state is then broadcasted by Peano to all other ranks.)
   */
  exahype::State _localState;

  /**
   * Flag indicating if one solver requested a local recomputation.
   * Is set in beginIteration(...).
   */
  static bool OneSolverRequestedLocalRecomputation;

  /**
   * A minimum time step size for each solver.
   */
  std::vector<double> _minTimeStepSizes;

  /**
   * A minimum cell size for each solver.
   */
  std::vector<double> _minCellSizes;

  /**
   * A maximum cell size for each solver.
   */
  std::vector<double> _maxCellSizes;

  /**
   * Prepare a appropriately sized vector _minTimeStepSizes
   * with elements initiliased to MAX_DOUBLE.
   */
  void prepareLocalTimeStepVariables();

  /**
   * Temporary variables for every registered
   * ADERDGSolver and LimitingADERDGSolver which are required for performing
   * the prediction.
   */
  exahype::solvers::PredictionTemporaryVariables _predictionTemporaryVariables;

  /**
   * Temporary variables for every registered solvers which are required
   * for performing the neighbour data merging.
   */
   exahype::solvers::MergingTemporaryVariables _mergingTemporaryVariables;

  /**
   * Initialises the temporary variables.
   *
   * \note We parallelise over the domain
   * (mapping is copied for each thread) and
   * over the solvers registered on a cell.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void initialiseTemporaryVariables();

  /**
   * Deletes the temporary variables.
   *
   * \note We need to initialise the temporary variables
   * in this mapping and not in the solvers since the
   * solvers in exahype::solvers::RegisteredSolvers
   * are not copied for every thread.
   */
  void deleteTemporaryVariables();

  /**
   * \return true if we perform a local recomputation for this solver.
   */
  static bool performLocalRecomputation(exahype::solvers::Solver* solver);

  /**
   * \return true if we have performed mesh refinement for this solver
   * and we do fused time stepping.
   *
   * \note That this is only necessary to check in \enterCell.
   * We will send out and reduce data for all solvers if we
   * run fused time stepping.
   */
  static bool performPrediction(exahype::solvers::Solver* solver);

  #ifdef Parallel
  /**
   * Drop incoming neighbour data for all solvers performing
   * a local recomputation.
   */
  static void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);

  /**
   * Merge incoming neighbour for all solvers performing
   * a local recomputation.
   */
  void mergeNeighourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS,int>&     src,
      const tarch::la::Vector<DIMENSIONS,int>&     dest,
      const int                                    destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);
  #endif

 public:
  /**
   * Run through the whole tree. Run concurrently on the fine grid.
   */
  peano::MappingSpecification enterCellSpecification(int level) const;
  /**
   * Run through the whole tree. Run concurrently on the fine grid.
   */
  peano::MappingSpecification leaveCellSpecification(int level) const;

  /**
   * Run through the whole tree. Avoid fine grid races.
   */
  peano::MappingSpecification touchVertexFirstTimeSpecification(int level) const;

  /**
   * Nop.
   */
  peano::MappingSpecification touchVertexLastTimeSpecification(int level) const;
  peano::MappingSpecification ascendSpecification(int level) const;
  peano::MappingSpecification descendSpecification(int level) const;
  peano::CommunicationSpecification communicationSpecification() const;

  /**
   * In case one solver requested a local recomputation,
   * advance limiter patches in time for this solvers.
   *
   * In case of fused time stepping, perform the prediction
   * for ADER-DG solvers which have performed mesh refinement.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * In case of fused time stepping, restrict data to next and top most (rank-local)
   * neighbour for solvers which have performed mesh refinement.
   *
   * \copydoc exahype::mappings::Prediction
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * In case one solver requested a local recomputation,
   * merge neighbouring limiter patches for this solvers.
   */
  void touchVertexFirstTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * Initialise temporary variables
   * if they are not initialised yet (or
   * if a new solver was introuced to the grid.
   * This is why we put the initialisation
   * in beginIteration().
   *
   * In debug mode, further resets counters for Riemann solves at
   * interior
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Frees previously allocated temporary variables.
   *
   * In debug mode, prints the output of counters.
   */
  void endIteration(exahype::State& solverState);

  #if defined(SharedMemoryParallelisation)
  /**
   * Copy the local state object over to the worker thread.
   */
  PredictionOrLocalRecomputation(const PredictionOrLocalRecomputation& masterThread);
  #endif

  /**
   * Free previously allocated temporary variables.
   */
  virtual ~PredictionOrLocalRecomputation();

#ifdef Parallel
  /**
   * In case one solver requested a local recomputation,
   * receive limiter boundary values from neighbouring ranks.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);
  /**
   * Send data out to the worker if fused time stepping is used.
   * We send for all solvers.
   * This does not depend on if we have performed a mesh
   * update or a local reecomputation.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);

  /**
   * Broadcast global data to a worker.
   *
   * \see exahype::Cell::broadcastGlobalDataToWorker
   */
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);

  /**
   * Receive global data to a from the master.
   *
   * \see exahype::Cell::mergeWithGlobalDataFromMaster
   */
  void receiveDataFromMaster(
      exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
      const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
      exahype::Vertex* const receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
      exahype::Cell& receivedCoarseGridCell,
      exahype::Vertex* const workersCoarseGridVertices,
      const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
      exahype::Cell& workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Reduce face data to the master for all solvers.
   * For solvers which requested a local recomputation,
   * reduce new time step data up to the master as well.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Merge with reduced face data for all solvers
   * from the workers.
   * For solvers which requested a local recomputation,
   * merge with new time step data from the workers.
   */
  void mergeWithMaster(
      const exahype::Cell& workerGridCell,
      exahype::Vertex* const workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker, const exahype::State& workerState,
      exahype::State& masterState);


  //
  // Below all methods are nop.
  //
  //===================================



  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);
  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);
  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);
  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);

  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);
  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif

#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  void mergeWithWorkerThread(const PredictionOrLocalRecomputation& workerThread);
#endif
  /**
   * Nop
   */
  PredictionOrLocalRecomputation();

  void createHangingVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * Nop.
   */
  void createInnerVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void createBoundaryVertex(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void destroyHangingVertex(
      const exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void destroyVertex(
      const exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void createCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void destroyCell(
      const exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  void touchVertexLastTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  /**
   * Nop.
   */
  void descend(
      exahype::Cell* const fineGridCells,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell);
  /**
   * Nop.
   */
  void ascend(exahype::Cell* const fineGridCells,
              exahype::Vertex* const fineGridVertices,
              const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
              exahype::Vertex* const coarseGridVertices,
              const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
              exahype::Cell& coarseGridCell);
};
#endif
