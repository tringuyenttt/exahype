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
 
#include "exahype/mappings/MeshRefinement.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "peano/grid/aspects/VertexStateAnalysis.h"

#include "tarch/la/VectorScalarOperations.h"

#include "tarch/multicore/Lock.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/LimiterStatusSpreading.h"

bool exahype::mappings::MeshRefinement::IsFirstIteration = true;

bool exahype::mappings::MeshRefinement::IsInitialMeshRefinement = false;

tarch::logging::Log exahype::mappings::MeshRefinement::_log("exahype::mappings::MeshRefinement");

tarch::multicore::BooleanSemaphore exahype::mappings::MeshRefinement::_semaphore;

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::MeshRefinement::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::MeshRefinement::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::MeshRefinement::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::MeshRefinement::MeshRefinement(const MeshRefinement& masterThread):
  _localState(masterThread._localState) {
}
#endif

void exahype::mappings::MeshRefinement::beginIteration(
  exahype::State& solverState
) {
  _localState = solverState;

  exahype::solvers::ADERDGSolver::Heap::getInstance().setName("ADERDGCellDescriptionHeap");
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().setName("FiniteVolumesCellDescriptionHeap");
  DataHeap::getInstance().setName("DataHeap");

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getMeshUpdateRequest()) {
      solver->zeroTimeStepSizes();
    }
    //assertion2(!solver->getNextMeshUpdateRequest(),solver->toString(),tarch::parallel::Node::getInstance().getRank());
  }

  #ifdef Parallel
  exahype::solvers::ADERDGSolver::Heap::getInstance().finishedToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();

  MetadataHeap::getInstance().setName("MetadataHeap");
  MetadataHeap::getInstance().finishedToSendSynchronousData();
  if (! MetadataHeap::getInstance().validateThatIncomingJoinBuffersAreEmpty() ) {
      exit(-1);
  }

  exahype::solvers::ADERDGSolver::Heap::getInstance().startToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().startToSendSynchronousData();
  DataHeap::getInstance().startToSendSynchronousData();
  MetadataHeap::getInstance().startToSendSynchronousData();
  #endif
}

void exahype::mappings::MeshRefinement::endIteration(exahype::State& solverState) {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getMeshUpdateRequest()) {
      solver->setNextAttainedStableState();
    }
  }

  #ifdef Parallel
  exahype::mappings::MeshRefinement::IsFirstIteration = false;
  #endif
}

void exahype::mappings::MeshRefinement::refineVertexIfNecessary(
  exahype::Vertex&                              fineGridVertex,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
  int                                           fineGridLevel,
  bool                                          isCalledByCreationalEvent
) const {
  bool refine =
      tarch::la::allGreater(
          fineGridH,exahype::solvers::Solver::getFinestMaximumMeshSizeOfAllSolvers());

  if (
    refine
    &&
    fineGridVertex.getRefinementControl()==Vertex::Records::Unrefined
  ) {
    switch ( _localState.mayRefine(isCalledByCreationalEvent,fineGridLevel) ) {
      case State::RefinementAnswer::DontRefineYet:
        break;
      case State::RefinementAnswer::Refine:
        fineGridVertex.refine();
        break;
      case State::RefinementAnswer::EnforceRefinement:
        fineGridVertex.enforceRefine();
        break;
    }
  }
}


void exahype::mappings::MeshRefinement::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  #ifdef Parallel
  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,false);
  #endif

  bool validCellDescriptionsIndexFound = false;
  dfor2(k)
    validCellDescriptionsIndexFound |= fineGridVertex.getCellDescriptionsIndex()[kScalar] > 0;
  enddforx

  if (validCellDescriptionsIndexFound) {
    dfor2(k)
       coarseGridVertices[ coarseGridVerticesEnumerator(k) ].vetoErasing();
    enddforx
  }

  bool hasToRefineToPreserveRegularityAlongBoundary = false;
  #ifdef Parallel
  bool boundaryVertexIsCritical =
    (
      (
      fineGridVertex.getRefinementControl()==Vertex::Records::Refining
      ||
      fineGridVertex.getRefinementControl()==Vertex::Records::RefinementTriggered
      )
      &&
      fineGridVertex.isBoundary()
      &&
      fineGridVertex.isRemote( _localState, true, false)
    );
  hasToRefineToPreserveRegularityAlongBoundary =
      coarseGridVerticesEnumerator.getLevel()+1 < peano::parallel::loadbalancing::Oracle::getInstance().getRegularLevelAlongBoundary() &&
      boundaryVertexIsCritical;
  #endif

  if (fineGridVertex.isMarkedForErasing() &&
      !hasToRefineToPreserveRegularityAlongBoundary &&
      fineGridVertex.getRefinementControl()==Vertex::Records::Refined
  ) {
    fineGridVertex.erase();
  }
}


void exahype::mappings::MeshRefinement::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createBoundaryVertex(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);

  logTraceOutWith1Argument("createBoundaryVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("createInnerVertex(...)", fineGridVertex, fineGridX,
                           fineGridH, coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  refineVertexIfNecessary(fineGridVertex,fineGridH,coarseGridVerticesEnumerator.getLevel()+1,true);

  logTraceOutWith1Argument("createInnerVertex(...)", fineGridVertex);
}


void exahype::mappings::MeshRefinement::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("createCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  fineGridCell.getCellData().setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);

  logTraceOutWith1Argument("createCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  fineGridVertex.mergeOnlyMetadata(_localState.getAlgorithmSection());

  bool erase =
      tarch::la::allSmallerEquals(
          fineGridH,exahype::solvers::Solver::getFinestMaximumMeshSizeOfAllSolvers());
  fineGridVertex.markForErasing(erase);

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::MeshRefinement::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  assertion(fineGridCell.isInside());

  bool refineFineGridCell = false;
  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if (solver->getMeshUpdateRequest()) {
      bool adjustSolution = exahype::mappings::MeshRefinement::IsFirstIteration;

      // Update limiter status and allocate limiter patch if necessary
      // Further evalute the limiter status based refinement criterion
      if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
        if (element!=exahype::solvers::Solver::NotFound) {
          const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>& indicesAdjacentToFineGridVertices =
              exahype::VertexOperations::readCellDescriptionsIndex(
                  fineGridVerticesEnumerator,fineGridVertices);
          if (multiscalelinkedcell::adjacencyInformationIsConsistent(
              indicesAdjacentToFineGridVertices)) {
            tarch::multicore::Lock lock(_semaphore);
            auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
            adjustSolution |=
                limitingADERDGSolver->updateLimiterStatusDuringLimiterStatusSpreading(
                    fineGridCell.getCellDescriptionsIndex(),element);
            refineFineGridCell |=
                limitingADERDGSolver->
                  evaluateLimiterStatusRefinementCriterion(
                        fineGridCell.getCellDescriptionsIndex(),element);
          }
        }
      }

      // TODO(Dominic): It is normally only necessary to check the
      // refinement criterion if we need to adjust the solution
      // However mark for refinement currently does some more things
      // Additionally, we would need to memorise where to refine /
      // erase in order to not break the behaviour of the
      // state machine
      refineFineGridCell |=
          solver->markForRefinement(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVertices,
              coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              IsInitialMeshRefinement,
              solverNumber);


      exahype::solvers::Solver::UpdateStateInEnterCellResult result =
          solver->updateStateInEnterCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVertices,
              coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              IsInitialMeshRefinement,
              solverNumber);
      refineFineGridCell |= result._refinementRequested;
      adjustSolution     |= result._newComputeCellAllocated;

      // Synchronise time stepping and adjust the solution
      // if necessary
      if (fineGridCell.isInitialised()) {
        const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
        if (element!=exahype::solvers::Solver::NotFound) {
          solver->zeroTimeStepSizes(fineGridCell.getCellDescriptionsIndex(),element);
          solver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),element);

          if (adjustSolution) {
            solver->adjustSolution(
                fineGridCell.getCellDescriptionsIndex(),
                element);

            if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
              static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                  updateLimiterStatusAndMinAndMaxAfterAdjustSolution(
                      fineGridCell.getCellDescriptionsIndex(),element);
            }
          }
        }

        exahype::Cell::resetNeighbourMergeFlags(
            fineGridCell.getCellDescriptionsIndex());
      }
    }
  }

  // Refine all adjacent vertices if necessary and possible.
  //    if (refineFineGridCell && _state.refineInitialGridInTouchVertexLastTime()) {
  if (refineFineGridCell) {
    dfor2(v)
      if (
          fineGridVertices[ fineGridVerticesEnumerator(v) ].getRefinementControl()==
          exahype::Vertex::Records::RefinementControl::Unrefined
          &&
          !fineGridVertices[ fineGridVerticesEnumerator(v) ].isHangingNode()
      ) {
        fineGridVertices[ fineGridVerticesEnumerator(v) ].refine();
      }
    enddforx
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::MeshRefinement::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);
  for (unsigned int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    // TODO(Dominic): Computationally heaviest operation is
    // restriction of volume unknowns up to a parent during
    // erasing operations. This is mostly localised.
    // Multicore-ising of this routine is currently postponed.

    if (solver->getMeshUpdateRequest()) {
          solver->updateStateInLeaveCell(
              fineGridCell,
              fineGridVertices,
              fineGridVerticesEnumerator,
              coarseGridCell,
              coarseGridVertices,
              coarseGridVerticesEnumerator,
              fineGridPositionOfCell,
              solverNumber);

      bool isStable = solver->attainedStableState(
          fineGridCell,
          fineGridVertices,
          fineGridVerticesEnumerator,
          solverNumber);

      solver->updateNextAttainedStableState(isStable);
    }
  }

  // shutdown metadata for empty cells (no cell descriptions)
  if (fineGridCell.isInitialised() &&
      fineGridCell.isEmpty()) {
    fineGridCell.shutdownMetaData();
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

#ifdef Parallel
void exahype::mappings::MeshRefinement::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  logTraceInWith6Arguments("mergeWithNeighbour(...)", vertex, neighbour,
                           fromRank, fineGridX, fineGridH, level);

  if (exahype::mappings::MeshRefinement::IsFirstIteration) {
    return;
  }
  vertex.mergeOnlyWithNeighbourMetadata(
      fromRank,fineGridX,fineGridH,level,
      _localState.getAlgorithmSection());

  logTraceOut("mergeWithNeighbour(...)");
}

void exahype::mappings::MeshRefinement::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  logTraceInWith5Arguments("prepareSendToNeighbour(...)", vertex,
                           toRank, x, h, level);

  vertex.sendOnlyMetadataToNeighbour(toRank,x,h,level);

  logTraceOut("prepareSendToNeighbour(...)");
}

bool exahype::mappings::MeshRefinement::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  if ( !_localState.isInvolvedInJoinOrFork() ) {
    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element =
          solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->prepareMasterCellDescriptionAtMasterWorkerBoundary(
            fineGridCell.getCellDescriptionsIndex(),element);
      }
    }
  }

  exahype::sendMasterWorkerCommunicationMetadata(
      worker,
      fineGridCell.getCellDescriptionsIndex(),
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getLevel());

  return true;
}

// TODO(Dominic): Add to docu: see documentation in peano/pdt/stdtemplates/MappingHeader.template
// on function receiveDataFromMaster
void exahype::mappings::MeshRefinement::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  receivedCell.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
  const int receivedMetadataIndex =
      exahype::receiveMasterWorkerCommunicationMetadata(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          receivedVerticesEnumerator.getCellCenter(),
          receivedVerticesEnumerator.getLevel());

  // Abusing the cell descriptions index
  receivedCell.setCellDescriptionsIndex(receivedMetadataIndex);
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if (localCell.isInitialised() &&
      !_localState.isInvolvedInJoinOrFork() ) {
    localCell.mergeWithMasterMetadata(
        receivedMasterCell.getCellDescriptionsIndex(),
        _localState.getAlgorithmSection());
  }

  MetadataHeap::getInstance().deleteData(receivedMasterCell.getCellDescriptionsIndex());
}

void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if (!localCell.isInside())
    return;

  if (localCell.isInitialised()) {
    exahype::solvers::ADERDGSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::sendCellDescriptions(toRank,localCell.getCellDescriptionsIndex(),
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
        auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      const int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if(element!=exahype::solvers::Solver::NotFound) {
        solver->sendDataToWorkerOrMasterDueToForkOrJoin(
            toRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
      } else {
        solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
      }
    }

    if (_localState.isJoiningWithMaster()) {
      exahype::solvers::ADERDGSolver::eraseCellDescriptions(localCell.getCellDescriptionsIndex());
      exahype::solvers::FiniteVolumesSolver::eraseCellDescriptions(localCell.getCellDescriptionsIndex());
      localCell.shutdownMetaData();
    }

  } else if (!localCell.isInitialised()){
    exahype::solvers::ADERDGSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);
    exahype::solvers::FiniteVolumesSolver::sendEmptyCellDescriptions(toRank,
        peano::heap::MessageType::ForkOrJoinCommunication,cellCentre,level);

    for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(toRank,cellCentre,level);
    }
  }
}

void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  if (!localCell.isInside())
    return;

  if (_localState.isNewWorkerDueToForkOfExistingDomain()) {
    localCell.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
    assertion1(!localCell.isInitialised(),localCell.toString());
  }
  if (!localCell.isInitialised()) {
    localCell.setupMetaData();
  }

  exahype::solvers::ADERDGSolver::mergeCellDescriptionsWithRemoteData(
      fromRank,localCell,
      peano::heap::MessageType::ForkOrJoinCommunication,
      cellCentre,level);
  exahype::solvers::FiniteVolumesSolver::mergeCellDescriptionsWithRemoteData(
      fromRank,localCell,
      peano::heap::MessageType::ForkOrJoinCommunication,
      cellCentre,level);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    const int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
          fromRank,localCell.getCellDescriptionsIndex(),element,cellCentre,level);
    } else {
      solver->dropWorkerOrMasterDataDueToForkOrJoin(fromRank,cellCentre,level);
    }
  }
}

void exahype::mappings::MeshRefinement::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    solver->sendMeshUpdateFlagsToMaster(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getLevel());

    if (!_localState.isInvolvedInJoinOrFork()) {
      const int element =
          solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->prepareWorkerCellDescriptionAtMasterWorkerBoundary(
            localCell.getCellDescriptionsIndex(),element);
      }
    }
  }

  exahype::sendMasterWorkerCommunicationMetadata(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      localCell.getCellDescriptionsIndex(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getLevel());
}

void exahype::mappings::MeshRefinement::mergeWithMaster(
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
  // Merge global solver states
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); solverNumber++) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    solver->mergeWithWorkerMeshUpdateFlags(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }

  // Merge cell states
  if (fineGridCell.isInitialised() &&
      _localState.isInvolvedInJoinOrFork()) {
    fineGridCell.mergeWithWorkerMetadata(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel(),
        _localState.getAlgorithmSection());
  } else {
    exahype::dropMetadata(
        worker,
        peano::heap::MessageType::MasterWorkerCommunication,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getLevel());
  }
}


//
// All methods below are nop,
//
// ==================================



void exahype::mappings::MeshRefinement::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::MeshRefinement::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::MeshRefinement::MeshRefinement() {
  // do nothing
}

exahype::mappings::MeshRefinement::~MeshRefinement() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
void exahype::mappings::MeshRefinement::mergeWithWorkerThread(
    const MeshRefinement& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::MeshRefinement::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::MeshRefinement::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::MeshRefinement::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::MeshRefinement::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
