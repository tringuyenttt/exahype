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
 
#include "exahype/mappings/TimeStepSizeComputation.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/parallel/Node.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include <limits>


bool exahype::mappings::TimeStepSizeComputation::SkipReductionInBatchedTimeSteps = false;



peano::CommunicationSpecification
exahype::mappings::TimeStepSizeComputation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterProcessingOfLocalSubtree,
      true);
}

peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
/**
 * Nop.
 */
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification exahype::mappings::TimeStepSizeComputation::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::TimeStepSizeComputation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::TimeStepSizeComputation::_log(
    "exahype::mappings::TimeStepSizeComputation");
int exahype::mappings::TimeStepSizeComputation::_mpiTag =
    tarch::parallel::Node::reserveFreeTag(
        "exahype::mappings::TimeStepSizeComputation");

void exahype::mappings::TimeStepSizeComputation::
    prepareEmptyLocalTimeStepData() {
  _minTimeStepSizes.resize(exahype::solvers::RegisteredSolvers.size());

  for (int i = 0;
       i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] = std::numeric_limits<double>::max();
  }
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation(
    const TimeStepSizeComputation& masterThread) {
  prepareEmptyLocalTimeStepData();
}

// Merge over threads
void exahype::mappings::TimeStepSizeComputation::mergeWithWorkerThread(
    const TimeStepSizeComputation& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], workerThread._minTimeStepSizes[i]);
  }
}
#endif

static double startNewTimeStepFV(
    const int cellDescriptionsIndex,
    const int element,
    const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices) {
  exahype::solvers::FiniteVolumesSolver::CellDescription& p =
      exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(cellDescriptionsIndex)[element];
  exahype::solvers::FiniteVolumesSolver* solver = static_cast<exahype::solvers::FiniteVolumesSolver*>(
      exahype::solvers::RegisteredSolvers[p.getSolverNumber()]);

  if (p.getType()==exahype::records::FiniteVolumesCellDescription::Cell) {
//         assertion1(p.getRefinementEvent()==exahype::records::FiniteVolumesCellDescription::None,p.toString()); // todo refine
    assertion1(multiscalelinkedcell::HangingVertexBookkeeper::allAdjacencyInformationIsAvailable(
        VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices)),fineGridVerticesEnumerator.toString());

    double* finiteVolumesSolutions[THREE_POWER_D];
    for (int nScalar=0; nScalar<THREE_POWER_D; ++nScalar) {
      if (exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(neighbourCellDescriptionsIndices[nScalar])) {
        exahype::records::FiniteVolumesCellDescription& pNeighbour =
            exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(neighbourCellDescriptionsIndices[nScalar])[p.getSolverNumber()]; // todo assumes same number of patches per cell
        finiteVolumesSolutions[nScalar] = exahype::DataHeap::getInstance().getData(pNeighbour.getSolution()).data();
      } else {
        finiteVolumesSolutions[nScalar] = exahype::DataHeap::getInstance().getData(p.getSolution()).data();
      }
    }

    double admissibleTimeStepSize = solver->stableTimeStepSize(
        finiteVolumesSolutions, p.getSize());

    assertion(!std::isnan(admissibleTimeStepSize));

    p.setTimeStamp(p.getTimeStamp()+p.getTimeStepSize());
    p.setTimeStepSize(admissibleTimeStepSize);

    return admissibleTimeStepSize;
  }

  return std::numeric_limits<double>::max();
}

void exahype::mappings::TimeStepSizeComputation::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    // ADER-DG
    const int numberOfSolvers = static_cast<int>(exahype::solvers::RegisteredSolvers.size());
    // please use a different UserDefined per mapping/event
    peano::datatraversal::autotuning::MethodTrace methodTrace = peano::datatraversal::autotuning::UserDefined2;
    int grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, methodTrace);

    pfor(solverNumber, 0, numberOfSolvers, grainSize)
      exahype::solvers::Solver* solver =
          exahype::solvers::RegisteredSolvers[solverNumber];
      int element = exahype::solvers::RegisteredSolvers[solverNumber]->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverNumber);

      if (element!=exahype::solvers::Solver::NotFound) {
        double admissibleTimeStepSize =
            solver->startNewTimeStep(fineGridCell.getCellDescriptionsIndex(),element);
        // TODO(Dominic): The FV solver does nothing at the moment in startNewTimeStep(...).
        // We currently rely on startNewTimeStepFV here. But this
        // function should be split in startNewTimeStep and updateSolution
        // solver functionality.
        if (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes) {
          const tarch::la::Vector<THREE_POWER_D, int> neighbourCellDescriptionsIndices = multiscalelinkedcell::getIndicesAroundCell(
                VertexOperations::readCellDescriptionsIndex(fineGridVerticesEnumerator, fineGridVertices));

          admissibleTimeStepSize =
              startNewTimeStepFV(fineGridCell.getCellDescriptionsIndex(),element,neighbourCellDescriptionsIndices);
        }
        _minTimeStepSizes[solverNumber] = std::min(
            admissibleTimeStepSize, _minTimeStepSizes[solverNumber]); // todo MPI
      }
    endpfor peano::datatraversal::autotuning::Oracle::getInstance()
        .parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::TimeStepSizeComputation::beginIteration(
    exahype::State& solverState) {
  prepareEmptyLocalTimeStepData();

  #ifdef Parallel
  DataHeap::getInstance().finishedToSendSynchronousData(); // See method documentation.
  DataHeap::getInstance().startToSendSynchronousData(); // See method documentation.
  #endif
}

void exahype::mappings::TimeStepSizeComputation::endIteration(
    exahype::State& solverState) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[i];

    logDebug("mergeLocalTimeStepDataIntoSolvers()",
             "solver " << i << " is updated with time step size "
                       << _minTimeStepSizes[i]);
    // TODO(Dominic): Fix bug:
    // This is buggy: The local solver already gets his
    // time step updated with an partially reduced time step size.
    // nextTimeStepSize must must be transferred to the solver before
    // we call startNewTimeStep() on the solver.
    solver->updateNextTimeStepSize(_minTimeStepSizes[i]);
    solver->startNewTimeStep();
  }
}

#ifdef Parallel
// Worker-master comm. Send to master.
void exahype::mappings::TimeStepSizeComputation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // TODO(Dominic): Merge with other Master-Worker exchanges: FaceUnknownsProjection

  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0,dt);
  }

  DataHeap::getInstance().sendData(
      _minTimeStepSizes.data(),_minTimeStepSizes.size(),
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getLevel(),
      peano::heap::MessageType::MasterWorkerCommunication);
}

// Worker-master comm. Receive from worker.
void exahype::mappings::TimeStepSizeComputation::mergeWithMaster(
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
    exahype::State& masterState) {
  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this loop if Debug/Asserts flag is not set.

  std::vector<double> receivedMinTimeStepSizes(_minTimeStepSizes.size());
  DataHeap::getInstance().receiveData(
      receivedMinTimeStepSizes.data(),receivedMinTimeStepSizes.size(),
      worker,
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getLevel(),
      peano::heap::MessageType::MasterWorkerCommunication);

  for (int i = 0; i < static_cast<int>(_minTimeStepSizes.size()); i++) {
    _minTimeStepSizes[i] =
        std::min(_minTimeStepSizes[i], receivedMinTimeStepSizes[i]);
  }

  for (auto dt : receivedMinTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this.
  for (auto dt : _minTimeStepSizes) {
    assertion1(dt>0.0,dt);
  } // Dead code elimination will get rid of this loop if Debug/Asserts flag is not set.
}


bool exahype::mappings::TimeStepSizeComputation::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  if (
    !peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
    &&
    SkipReductionInBatchedTimeSteps
  ) {
    return false;
  }
  else return true;
}


void exahype::mappings::TimeStepSizeComputation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::TimeStepSizeComputation::TimeStepSizeComputation() {
  // do nothing
}

exahype::mappings::TimeStepSizeComputation::~TimeStepSizeComputation() {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::TimeStepSizeComputation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
