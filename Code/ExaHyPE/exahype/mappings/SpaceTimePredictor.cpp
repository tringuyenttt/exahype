#include "exahype/mappings/SpaceTimePredictor.h"

#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "peano/datatraversal/autotuning/Oracle.h"

#include "exahype/solvers/Solver.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::SpaceTimePredictor::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}
peano::MappingSpecification
exahype::mappings::SpaceTimePredictor::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log exahype::mappings::SpaceTimePredictor::_log( "exahype::mappings::SpaceTimePredictor");
int                 exahype::mappings::SpaceTimePredictor::_mpiTag = tarch::parallel::Node::reserveFreeTag( "exahype::mappings::SpaceTimePredictor" );



exahype::mappings::SpaceTimePredictor::SpaceTimePredictor() {
  // do nothing
}

exahype::mappings::SpaceTimePredictor::~SpaceTimePredictor() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SpaceTimePredictor::SpaceTimePredictor(
  const SpaceTimePredictor& masterThread
) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorkerThread(
    const SpaceTimePredictor& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::SpaceTimePredictor::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::destroyCell(
    const exahype::Cell& fineGridCell,
exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::SpaceTimePredictor::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::prepareSendToNeighbour(
    exahype::Vertex&                              vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const tarch::la::Vector<DIMENSIONS, double>&  h,
    int                                           level
) {
  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getADERDGCellDescriptionsIndex();

  /* Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
   * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
   * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
   *
   * Note that from the viewpoint of a cell, the face
   * has always the "opposite" index, i.e., we solve a Riemann
   * problem on the left face of the right cell (which
   * is the right face of the left cell).
   */
  constexpr int cellIndicesLeft  [4] = {0, 2, 4, 6};
  constexpr int cellIndicesRight [4] = {1, 3, 5, 7};
  constexpr int cellIndicesFront [4] = {0, 1, 4, 5};
  constexpr int cellIndicesBack  [4] = {2, 3, 6, 7};
#if DIMENSIONS == 3
  constexpr int cellIndicesBottom[4] = {0, 1, 2, 3};
  constexpr int cellIndicesTop   [4] = {4, 5, 6, 7};
#endif
  dfor2(dest)
  dfor2(src)
    if (
      vertex.getAdjacentRanks()(destScalar)==toRank
      &&
      vertex.getAdjacentRanks()(srcScalar)==tarch::parallel::Node::getInstance().getRank()
      &&
      tarch::la::countEqualEntries(dest,src)==1    // we are solely exchanging faces
    ) {
      std::vector<records::ADERDGCellDescription>& cellDescriptions = ADERDGCellDescriptionHeap::getInstance().getData(srcScalar);
      for (int currentSolver=0; currentSolver<static_cast<int>(cellDescriptions.size()); currentSolver++) {
        if (cellDescriptions[currentSolver].getType()==exahype::records::ADERDGCellDescription::RealCell) {
          exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescriptions[currentSolver].getSolverNumber()];

          const int numberOfFaceDof       = solver->getUnknownsPerFace();
          const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src,dest);
          assertion(normalOfExchangedFace>=0 && normalOfExchangedFace<DIMENSIONS);
          const int offsetInFaceArray     = 2*normalOfExchangedFace + src(normalOfExchangedFace)<dest(normalOfExchangedFace) ? 1 : 0;

          const double* lQhbnd = DataHeap::getInstance().getData(cellDescriptions[currentSolver].getExtrapolatedPredictor()).data() + (offsetInFaceArray * numberOfFaceDof);
          const double* lFhbnd = DataHeap::getInstance().getData(cellDescriptions[currentSolver].getFluctuation()).data()           + (offsetInFaceArray * numberOfFaceDof);

          DataHeap::getInstance().sendData( lQhbnd, numberOfFaceDof, toRank, x, level, peano::heap::MessageType::NeighbourCommunication );
          DataHeap::getInstance().sendData( lFhbnd, numberOfFaceDof, toRank, x, level, peano::heap::MessageType::NeighbourCommunication );
        }
        else {
          assertionMsg( false, "Dominic, please implement" );
        }
      }
    }
  enddforx
  enddforx
}


void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}


bool exahype::mappings::SpaceTimePredictor::prepareSendToWorker(
  exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
  const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
  exahype::Vertex* const coarseGridVertices,
  const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
  exahype::Cell&                            coarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
  int worker
) {
  for (
    std::vector<exahype::solvers::Solver*>::iterator p = exahype::solvers::RegisteredSolvers.begin();
    p != exahype::solvers::RegisteredSolvers.end();
    p++
  ) {
    (*p)->sendToRank(worker, _mpiTag);
  }

  return true;
}


void exahype::mappings::SpaceTimePredictor::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithMaster(
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
  // do nothing
}


void exahype::mappings::SpaceTimePredictor::receiveDataFromMaster(
  exahype::Cell&                              receivedCell,
  exahype::Vertex*                            receivedVertices,
  const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
  exahype::Vertex* const                      receivedCoarseGridVertices,
  const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
  exahype::Cell&                              receivedCoarseGridCell,
  exahype::Vertex* const                      workersCoarseGridVertices,
  const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
  exahype::Cell&                              workersCoarseGridCell,
  const tarch::la::Vector<DIMENSIONS, int>&   fineGridPositionOfCell
) {
  for (
    std::vector<exahype::solvers::Solver*>::iterator p = exahype::solvers::RegisteredSolvers.begin();
    p != exahype::solvers::RegisteredSolvers.end();
    p++
  ) {
    (*p)->receiveFromRank(tarch::parallel::NodePool::getInstance().getMasterRank(),_mpiTag);
  }
}


void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::SpaceTimePredictor::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (ADERDGCellDescriptionHeap::getInstance().
      isValidIndex(fineGridCell.getADERDGCellDescriptionsIndex())) {
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance()
        .getData(fineGridCell.getADERDGCellDescriptionsIndex())
        .size());

    // please use a different UserDefined per mapping/event
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined7;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);
    pfor(i, 0, numberOfADERDGCellDescriptions, grainSize)
      records::ADERDGCellDescription& p =
          fineGridCell.getADERDGCellDescription(i);

      if (p.getType()==exahype::records::ADERDGCellDescription::RealCell) {
        exahype::solvers::Solver* solver =
            exahype::solvers::RegisteredSolvers[p.getSolverNumber()];

        // space-time DoF (basisSize**(DIMENSIONS+1))
        double* lQi =
            DataHeap::getInstance().
            getData(p.getSpaceTimePredictor()).
            data();
        double* lFi =
            DataHeap::getInstance().
            getData(p.getSpaceTimeVolumeFlux()).
            data();

        // volume DoF (basisSize**(DIMENSIONS))
        double* luh = DataHeap::getInstance().
            getData(p.getSolution()).
            data();
        double* lQhi = DataHeap::getInstance().
            getData(p.getPredictor()).
            data();
        double* lFhi = DataHeap::getInstance().
            getData(p.getVolumeFlux()).
            data();

        // face DoF (basisSize**(DIMENSIONS-1))
        double* lQhbnd = DataHeap::getInstance().
            getData(p.getExtrapolatedPredictor()).
            data();
        double* lFhbnd = DataHeap::getInstance().
            getData(p.getFluctuation()).
            data();

        solver->spaceTimePredictor(
            lQi, lFi, lQhi, lFhi,
            lQhbnd,  // da kommt was drauf todo what does this mean?
            lFhbnd,  // da kommt was drauf todo what does this mean?
            luh, fineGridVerticesEnumerator.getCellSize(),
            p.getPredictorTimeStepSize());

        logDebug("enterCell(...)", "debug::after::luh[0]" << luh[0]);
        logDebug("enterCell(...)", "debug::after::lQi[0]" << lQi[0]);
        logDebug("enterCell(...)", "debug::after::lFi[0]" << lFi[0]);
        logDebug("enterCell(...)", "debug::after::lQhi[0]" << lQhi[0]);
        logDebug("enterCell(...)", "debug::after::lFhi[0]" << lFhi[0]);
        logDebug("enterCell(...)", "debug::after::lQhbnd[0]" << lQhbnd[0]);
        logDebug("enterCell(...)", "debug::after::lFhbnd[0]" << lFhbnd[0]);
      }
    endpfor

    peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SpaceTimePredictor::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::beginIteration(exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SpaceTimePredictor::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
