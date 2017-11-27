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
 
#include "exahype/mappings/LocalRecomputationOrPrediction.h"

#include <algorithm>

#include "tarch/multicore/Loop.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

tarch::logging::Log exahype::mappings::LocalRecomputationOrPrediction::_log(
    "exahype::mappings::LocalRecomputationOrPrediction");

peano::CommunicationSpecification
exahype::mappings::LocalRecomputationOrPrediction::communicationSpecification() const {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::enterCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::leaveCellSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::touchVertexFirstTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

// Below specs are all nop
peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::ascendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::LocalRecomputationOrPrediction::descendSpecification(int level) const {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}


void exahype::mappings::LocalRecomputationOrPrediction::prepareLocalTimeStepVariables(){
  const unsigned int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  _minTimeStepSizes.resize(numberOfSolvers);
  _minCellSizes.resize(numberOfSolvers);
  _maxCellSizes.resize(numberOfSolvers);

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] = std::numeric_limits<double>::max();
    _minCellSizes    [solverNumber] = std::numeric_limits<double>::max();
    _maxCellSizes    [solverNumber] = -std::numeric_limits<double>::max(); // "-", min
  }
}

void exahype::mappings::LocalRecomputationOrPrediction::initialiseTemporaryVariables() {
  exahype::solvers::initialiseTemporaryVariables(_predictionTemporaryVariables);
  exahype::solvers::initialiseTemporaryVariables(_mergingTemporaryVariables);
}

void exahype::mappings::LocalRecomputationOrPrediction::deleteTemporaryVariables() {
  exahype::solvers::deleteTemporaryVariables(_predictionTemporaryVariables);
  exahype::solvers::deleteTemporaryVariables(_mergingTemporaryVariables);
}


exahype::mappings::LocalRecomputationOrPrediction::LocalRecomputationOrPrediction()
  #ifdef Debug
  :
  _interiorFaceMerges(0),
  _boundaryFaceMerges(0)
  #endif
{
  // do nothing
}

exahype::mappings::LocalRecomputationOrPrediction::~LocalRecomputationOrPrediction() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::LocalRecomputationOrPrediction::LocalRecomputationOrPrediction(
    const LocalRecomputationOrPrediction& masterThread)
: _localState(masterThread._localState) {
  prepareLocalTimeStepVariables();

  initialiseTemporaryVariables();
}
// Merge over threads
void exahype::mappings::LocalRecomputationOrPrediction::mergeWithWorkerThread(
    const LocalRecomputationOrPrediction& workerThread) {
  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _minTimeStepSizes[solverNumber] =
        std::min(_minTimeStepSizes[solverNumber], workerThread._minTimeStepSizes[solverNumber]);
    _minCellSizes[solverNumber] =
        std::min(_minCellSizes[solverNumber], workerThread._minCellSizes[solverNumber]);
    _maxCellSizes[solverNumber] =
        std::max(_maxCellSizes[solverNumber], workerThread._maxCellSizes[solverNumber]);
  }
}
#endif

void exahype::mappings::LocalRecomputationOrPrediction::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  prepareLocalTimeStepVariables();

  initialiseTemporaryVariables();

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceMerges = 0;
  _boundaryFaceMerges = 0;
  #endif

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

bool exahype::mappings::LocalRecomputationOrPrediction::performLocalRecomputation(
    exahype::solvers::Solver* solver) {
  return
      solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
      &&
      static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainChange()
      ==exahype::solvers::LimiterDomainChange::Irregular;
}

bool exahype::mappings::LocalRecomputationOrPrediction::performPrediction(
    exahype::solvers::Solver* solver) {
  return exahype::State::fuseADERDGPhases() &&
         solver->getMeshUpdateRequest();
}

void exahype::mappings::LocalRecomputationOrPrediction::endIteration(
    exahype::State& state) {
  logTraceInWith1Argument("endIteration(State)", state);

  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
    if ( performLocalRecomputation( solver ) ) {
      logDebug("endIteration(state)","_minCellSizes[solverNumber]="<<_minCellSizes[solverNumber]<<
          ",_minCellSizes[solverNumber]="<<_maxCellSizes[solverNumber]);
      assertion1(std::isfinite(_minTimeStepSizes[solverNumber]),_minTimeStepSizes[solverNumber]);
      assertion1(_minTimeStepSizes[solverNumber]>0.0,_minTimeStepSizes[solverNumber]);

      solver->updateNextMinCellSize(_minCellSizes[solverNumber]);
      solver->updateNextMaxCellSize(_maxCellSizes[solverNumber]);
      if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
        assertion3(solver->getNextMinCellSize()<std::numeric_limits<double>::max(),
            solver->getNextMinCellSize(),_minCellSizes[solverNumber],solver->toString());
        assertion3(solver->getNextMaxCellSize()>0,
            solver->getNextMaxCellSize(),_maxCellSizes[solverNumber],solver->toString());
      }

      solver->updateMinNextTimeStepSize(_minTimeStepSizes[solverNumber]);

      if (
          exahype::State::fuseADERDGPhases()
          #ifdef Parallel
          && tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getGlobalMasterRank()
          #endif
      ) {
        exahype::solvers::Solver::
        reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(solver);
      }
      if (exahype::State::fuseADERDGPhases()) {
        solver->startNewTimeStepFused(true,true);
      } else {
        solver->startNewTimeStep();
      }

      logDebug("endIteration(state)","updatedTimeStepSize="<<solver->getMinTimeStepSize());
    }
  }

  deleteTemporaryVariables();

  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interiorFaceSolves: " << _interiorFaceMerges);
  logInfo("endIteration(...)","boundaryFaceSolves: " << _boundaryFaceMerges);
  #endif

  logTraceOutWith1Argument("endIteration(State)", state);
}

void exahype::mappings::LocalRecomputationOrPrediction::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined3);
    pfor(solverNumber, 0, numberOfSolvers, grainSize.getGrainSize())
      auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if ( element!=exahype::solvers::Solver::NotFound ) {
        if ( performLocalRecomputation( solver ) ) {
          auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
          limitingADERDG->recomputeSolutionLocally(
              fineGridCell.getCellDescriptionsIndex(),element);

          double admissibleTimeStepSize = std::numeric_limits<double>::max();
          if (exahype::State::fuseADERDGPhases()) {
            limitingADERDG->recomputePredictorLocally(
                fineGridCell.getCellDescriptionsIndex(),element,
                _predictionTemporaryVariables);
            admissibleTimeStepSize = limitingADERDG->startNewTimeStepFused(
                fineGridCell.getCellDescriptionsIndex(),element,
                exahype::State::isFirstIterationOfBatchOrNoBatch(),
                exahype::State::isLastIterationOfBatchOrNoBatch());
          } else {
            admissibleTimeStepSize = limitingADERDG->startNewTimeStep(
                fineGridCell.getCellDescriptionsIndex(),element);
          }

          _minTimeStepSizes[solverNumber] = std::min(
              admissibleTimeStepSize, _minTimeStepSizes[solverNumber]);
          _minCellSizes[solverNumber] = std::min(
              fineGridVerticesEnumerator.getCellSize()[0],_minCellSizes[solverNumber]);
          _maxCellSizes[solverNumber] = std::max(
              fineGridVerticesEnumerator.getCellSize()[0],_maxCellSizes[solverNumber]);

          limitingADERDG->determineMinAndMax(fineGridCell.getCellDescriptionsIndex(),element);
        } else if (
            exahype::State::fuseADERDGPhases() &&
            solver->getMeshUpdateRequest() ) {
          exahype::solvers::ADERDGSolver::performPredictionAndVolumeIntegral(
              solver,fineGridCell.getCellDescriptionsIndex(),element,_predictionTemporaryVariables);

          solver->prolongateDataAndPrepareDataRestriction(fineGridCell.getCellDescriptionsIndex(),element);
        }
      }
    endpfor
    grainSize.parallelSectionHasTerminated();

    exahype::Cell::resetNeighbourMergeFlags(
        fineGridCell.getCellDescriptionsIndex());
    exahype::Cell::resetFaceDataExchangeCounters(
        fineGridCell.getCellDescriptionsIndex(),
        fineGridVertices,fineGridVerticesEnumerator);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}


void exahype::mappings::LocalRecomputationOrPrediction::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  exahype::mappings::Prediction::restrictDataAndPostProcess(
      fineGridCell,coarseGridCell,
      exahype::State::Records::AlgorithmSection::MeshRefinementOrLocalOrGlobalRecomputationAllSend);

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}

void exahype::mappings::LocalRecomputationOrPrediction::touchVertexFirstTime(
  exahype::Vertex& fineGridVertex,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
  const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
  exahype::Cell& coarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  dfor2(pos1)
    dfor2(pos2)
      if (fineGridVertex.hasToMergeNeighbours(pos1,pos1Scalar,pos2,pos2Scalar,fineGridX,fineGridH)) { // Assumes that we have to valid indices
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined5);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
          auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
          if ( performLocalRecomputation(solver) ) {
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
                  mergeNeighboursBasedOnLimiterStatus(
                  cellDescriptionsIndex1,element1,
                  cellDescriptionsIndex2,element2,
                  pos1,pos2,
                  true, /* isRecomputation */
                  _mergingTemporaryVariables._tempFaceUnknowns[solverNumber]);
            }
          }

          #ifdef Debug // TODO(Dominic)
          _interiorFaceMerges++;
          #endif
        endpfor
        grainSize.parallelSectionHasTerminated();

        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
      if (fineGridVertex.hasToMergeWithBoundaryData(pos1,pos1Scalar,pos2,pos2Scalar,fineGridX,fineGridH)) {
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined6);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
          auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];
          const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
          const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
          int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
          int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
          assertion4((element1==exahype::solvers::Solver::NotFound &&
                      element2==exahype::solvers::Solver::NotFound)
                     || (element1 >= 0 && element2==exahype::solvers::Solver::NotFound)
                     || (element2 >= 0 && element1==exahype::solvers::Solver::NotFound),
                     cellDescriptionsIndex1,cellDescriptionsIndex2,element1,element2); // TODO(Dominic): Move down
          if (
              element1 >= 0 &&
              performLocalRecomputation(solver)
          ) {
            auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
            auto& solverPatch1 = limitingADERDG->getSolver()->getCellDescription(cellDescriptionsIndex1,element1);

            limitingADERDG->
              mergeWithBoundaryDataBasedOnLimiterStatus( // !!! Be aware of indices "2" and "1" and the order of the arguments.
                cellDescriptionsIndex1,element1,
                solverPatch1.getLimiterStatus(), // !!! We assume here that we have already unified the merged limiter status values.
                pos1,pos2,                              // The cell-based limiter status is still holding the old value though.
                true,
                _mergingTemporaryVariables._tempFaceUnknowns[solverNumber]);

            #ifdef Debug
            _boundaryFaceMerges++;
            #endif
          }
          else if (
              element2 >= 0 &&
              performLocalRecomputation(solver)
          ){
            auto* limitingADERDG = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
            auto& solverPatch2 = limitingADERDG->getSolver()->getCellDescription(cellDescriptionsIndex2,element2);

            limitingADERDG->
              mergeWithBoundaryDataBasedOnLimiterStatus( // !!! Be aware of indices "2" and "1" and the order of the arguments.
                cellDescriptionsIndex2,element2,
                solverPatch2.getLimiterStatus(), // !!! We assume here that we have already unified the merged limiter status values
                pos2,pos1,                              // The cell-based limiter status is still holding the old value though.
                true,
                _mergingTemporaryVariables._tempFaceUnknowns[solverNumber]);
            #ifdef Debug
            _boundaryFaceMerges++;
            #endif
          }
        endpfor
        grainSize.parallelSectionHasTerminated();

        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
    enddforx
  enddforx
}


#ifdef Parallel
///////////////////////////////////////
// NEIGHBOUR
///////////////////////////////////////
void exahype::mappings::LocalRecomputationOrPrediction::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  if (!vertex.hasToCommunicate(fineGridH)) {
    return;
  }

  dfor2(myDest)
    dfor2(mySrc)
      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

      int destScalar = TWO_POWER_D - myDestScalar - 1;
      int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

      if (vertex.hasToReceiveMetadata(fromRank,src,dest)) {
        const int receivedMetadataIndex =
            exahype::receiveNeighbourCommunicationMetadata(
                fromRank, fineGridX, level);
        exahype::MetadataHeap::HeapEntries& receivedMetadata =
            MetadataHeap::getInstance().getData(receivedMetadataIndex);
        assertionEquals(receivedMetadata.size(),
            exahype::NeighbourCommunicationMetadataPerSolver*solvers::RegisteredSolvers.size());

        if(vertex.hasToMergeWithNeighbourData(src,dest)) { // Only communicate data once per face
          mergeNeighourData(
              fromRank,
              src,dest,
              vertex.getCellDescriptionsIndex()[destScalar],
              fineGridX,level,
              receivedMetadata);

          vertex.setFaceDataExchangeCountersOfDestination(src,dest,TWO_POWER_D); // !!! Do not forget this
          vertex.setMergePerformed(src,dest,true);
        } else {
          dropNeighbourData(
              fromRank,
              src,dest,
              fineGridX,level,
              receivedMetadata);
        }
        // Clean up
        MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
      }
    enddforx
  enddforx
}

void exahype::mappings::LocalRecomputationOrPrediction::dropNeighbourData(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];

    if ( performLocalRecomputation( solver ) ) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

      logDebug("dropNeighbourData(...)", "drop data for solver " << solverNumber << " from rank " <<
              fromRank << " at vertex x=" << x << ", level=" << level <<
              ", src=" << src << ", dest=" << dest);

      limitingADERDGSolver->dropNeighbourSolverAndLimiterData(fromRank,src,dest,x,level);
    }
  }
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeNeighourData(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level,
    const exahype::MetadataHeap::HeapEntries&    receivedMetadata) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(destCellDescriptionIndex));

  for(unsigned int solverNumber = solvers::RegisteredSolvers.size(); solverNumber-- > 0;) {
    auto* solver = solvers::RegisteredSolvers[solverNumber];

    if ( performLocalRecomputation( solver ) ) {
      const int element = solver->tryGetElement(destCellDescriptionIndex,solverNumber);
      const int offset  = exahype::NeighbourCommunicationMetadataPerSolver*solverNumber;

      if (element!=exahype::solvers::Solver::NotFound
          && receivedMetadata[offset].getU()!=exahype::InvalidMetadataEntry) {
        logDebug("mergeWithNeighbourData(...)", "receive data for solver " << solverNumber << " from rank " <<
                      fromRank << " at vertex x=" << x << ", level=" << level <<
                      ", src=" << src << ", dest=" << dest);

        exahype::MetadataHeap::HeapEntries metadataPortion(
                  receivedMetadata.begin()+offset,
                  receivedMetadata.begin()+offset+exahype::NeighbourCommunicationMetadataPerSolver);

        auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDGSolver->mergeWithNeighbourDataBasedOnLimiterStatus(
            fromRank,
            metadataPortion,
            destCellDescriptionIndex,element,src,dest,
            true, /* isRecomputation */
            _mergingTemporaryVariables._tempFaceUnknowns[solverNumber],
            x,level);
      } else {
        logDebug("mergeWithNeighbourData(...)", "drop data for solver " << solverNumber << " from rank " <<
                      fromRank << " at vertex x=" << x << ", level=" << level <<
                      ", src=" << src << ", dest=" << dest);

        auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
        limitingADERDGSolver->dropNeighbourSolverAndLimiterData(fromRank,src,dest,x,level);
      }
    }
  }
}

void exahype::mappings::LocalRecomputationOrPrediction::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  if ( exahype::State::fuseADERDGPhases() ) {
    vertex.sendToNeighbour(toRank,x,h,level);
  }
}

bool exahype::mappings::LocalRecomputationOrPrediction::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  exahype::Cell::broadcastGlobalDataToWorker(
      worker,
      coarseGridVerticesEnumerator.getCellCenter(),
      coarseGridVerticesEnumerator.getLevel());

  return true;
}

void exahype::mappings::LocalRecomputationOrPrediction::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  exahype::Cell::mergeWithGlobalDataFromMaster(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      workersCoarseGridVerticesEnumerator.getCellCenter(),
      workersCoarseGridVerticesEnumerator.getLevel());
}

void exahype::mappings::LocalRecomputationOrPrediction::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  exahype::Cell::reduceGlobalDataToMaster(
      tarch::parallel::NodePool::getInstance().getMasterRank(),
      verticesEnumerator.getCellCenter(),
      verticesEnumerator.getLevel());

  if ( exahype::State::fuseADERDGPhases() ) {
    localCell.reduceDataToMasterPerCell(
        tarch::parallel::NodePool::getInstance().getMasterRank(),
        verticesEnumerator.getCellCenter(),
        verticesEnumerator.getCellSize(),
        verticesEnumerator.getLevel());
  }
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeWithMaster(
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
  exahype::Cell::mergeWithGlobalDataFromWorker(
      worker,
      fineGridVerticesEnumerator.getCellCenter(),
      fineGridVerticesEnumerator.getLevel());

  if ( exahype::State::fuseADERDGPhases() ) {
    fineGridCell.mergeWithDataFromWorkerPerCell(
        worker,
        fineGridVerticesEnumerator.getCellCenter(),
        fineGridVerticesEnumerator.getCellSize(),
        fineGridVerticesEnumerator.getLevel());
  }
}

//
// Below all methods are nop.
//
//=====================================


void exahype::mappings::LocalRecomputationOrPrediction::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::LocalRecomputationOrPrediction::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::LocalRecomputationOrPrediction::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
