#include "exahype/mappings/RiemannSolver.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Lock.h"
#include "tarch/multicore/Loop.h"

#include "exahype/solvers/Solver.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::CommunicationSpecification
exahype::mappings::RiemannSolver::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          SendDataAndStateAfterLastTouchVertexLastTime,
      true);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::RiemannSolver::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces);
}

tarch::logging::Log exahype::mappings::RiemannSolver::_log(
    "exahype::mappings::RiemannSolver");

exahype::mappings::RiemannSolver::RiemannSolver() {
  // do nothing
}

exahype::mappings::RiemannSolver::~RiemannSolver() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::RiemannSolver::RiemannSolver(
    const RiemannSolver& masterThread)
    : _localState(masterThread._localState) {}

void exahype::mappings::RiemannSolver::mergeWithWorkerThread(
    const RiemannSolver& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::RiemannSolver::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::RiemannSolver::mergeWithNeighbour(
  exahype::Vertex&                              vertex,
  const exahype::Vertex&                        neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridX,
  const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
  int                                           level) {
  #if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
  #endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getADERDGCellDescriptionsIndex();

  dfor2(dest)
    dfor2(src)
      if (
        vertex.getAdjacentRanks()(destScalar) ==
        tarch::parallel::Node::getInstance().getRank() &&
        vertex.getAdjacentRanks()(srcScalar) == fromRank &&
        tarch::la::countEqualEntries(dest, src) == 1  // we are solely exchanging faces
      ) {
        const int destCellDescriptionIndex =
          adjacentADERDGCellDescriptionsIndices(destScalar);

        assertion5(
          ADERDGCellDescriptionHeap::getInstance().isValidIndex(
          destCellDescriptionIndex),
          src, dest, multiscalelinkedcell::indicesToString(adjacentADERDGCellDescriptionsIndices),
          vertex.toString(),
          tarch::parallel::Node::getInstance().getRank()
        );

    std::vector<records::ADERDGCellDescription>& cellDescriptions =
        ADERDGCellDescriptionHeap::getInstance().getData(
            destCellDescriptionIndex);

    for (
      int currentSolver = 0;
      currentSolver < static_cast<int>(cellDescriptions.size());
      currentSolver++
    ) {
      if (cellDescriptions[currentSolver].getType() == exahype::records::ADERDGCellDescription::Cell) {
        exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescriptions[currentSolver].getSolverNumber()];

        const int numberOfFaceDof = solver->getUnknownsPerFace();
        const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
        assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
        const int offsetInFaceArray =
            2 * normalOfExchangedFace + src(normalOfExchangedFace) < dest(normalOfExchangedFace) ?
            0 : 1;

        assertion(DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getExtrapolatedPredictor()));
        assertion(DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getFluctuation()));

        const double* lQhbnd =
            DataHeap::getInstance().getData(cellDescriptions[currentSolver].getExtrapolatedPredictor()).data() +
            (offsetInFaceArray * numberOfFaceDof);
        const double* lFhbnd =
            DataHeap::getInstance().getData(cellDescriptions[currentSolver].getFluctuation()).data() +
            (offsetInFaceArray * numberOfFaceDof);

        if (adjacentADERDGCellDescriptionsIndices(destScalar) == multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex) {
          #if defined(PeriodicBC)
          assertionMsg(false, "Vasco, we have to implement this");
          #else
          assertionMsg(false, "should never been entered");
          #endif
        } else {
          logDebug(
            "mergeWithNeighbour(...)",
            "receive two arrays from rank " << fromRank << " for vertex " << vertex.toString()
            << ", src type=" << multiscalelinkedcell::indexToString(adjacentADERDGCellDescriptionsIndices(srcScalar))
          );

          int receivedlQhbndIndex = DataHeap::getInstance().createData(0,numberOfFaceDof);
          int receivedlFhbndIndex = DataHeap::getInstance().createData(0,numberOfFaceDof);

          assertion( DataHeap::getInstance().getData(receivedlQhbndIndex).empty() );
          assertion( DataHeap::getInstance().getData(receivedlFhbndIndex).empty() );

          DataHeap::getInstance().receiveData(
              receivedlQhbndIndex,
              fromRank, fineGridX, level,
              peano::heap::MessageType::NeighbourCommunication);
          DataHeap::getInstance().receiveData(
              receivedlFhbndIndex,
              fromRank, fineGridX, level,
              peano::heap::MessageType::NeighbourCommunication);

          switch (normalOfExchangedFace) {
            case 0:
              break;
            case 1:
              break;
            case 2:
              break;
            default:
              assertionMsg(false, "we do not support 4d");
              break;
          }
/*
          solveRiemannProblemAtInterface(
            adjacentCellDescriptionsIndices[cellIndicesLeft[i]],
            adjacentCellDescriptionsIndices[cellIndicesRight[i]],
            EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT, 0);

          solveRiemannProblemAtInterface(
            adjacentCellDescriptionsIndices[cellIndicesFront[i]],
            adjacentCellDescriptionsIndices[cellIndicesBack[i]],
            EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT, 0);

          #if DIMENSIONS == 3
          solveRiemannProblemAtInterface(
            adjacentCellDescriptionsIndices[cellIndicesBottom[i]],
            adjacentCellDescriptionsIndices[cellIndicesTop[i]],
            EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM, 0);
          #endif
*/

          DataHeap::getInstance().deleteData(receivedlQhbndIndex);
          DataHeap::getInstance().deleteData(receivedlFhbndIndex);
        }
      } else {
        assertionMsg(false, "Dominic, please implement");
      }
    }
  }
  enddforx enddforx
}

void exahype::mappings::RiemannSolver::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::RiemannSolver::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing but please consult header documentation.

  return true;
}

void exahype::mappings::RiemannSolver::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithMaster(
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

void exahype::mappings::RiemannSolver::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing but please consult header documentation
}

void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::RiemannSolver::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::RiemannSolver::touchVertexFirstTime(
    exahype::Vertex&                              fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>&  fineGridH,
    exahype::Vertex* const                        coarseGridVertices,
    const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
    exahype::Cell&                                coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>&     fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                           fineGridX, fineGridH,
                           coarseGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfVertex);

  tarch::la::Vector<TWO_POWER_D, int>& adjacentCellDescriptionsIndices =
      fineGridVertex.getADERDGCellDescriptionsIndex();
  logDebug(
      "touchVertexFirstTime(...)",
      "cell descriptions around vertex. "
          << "coarse grid level: " << coarseGridVerticesEnumerator.getLevel()
          << ", fine grid position:" << fineGridPositionOfVertex
          << ", adjacent cell descriptions indices:"
          << adjacentCellDescriptionsIndices);
  logDebug("touchVertexFirstTime(...)", "cell descriptions around vertex. "
                                            << "fine grid x " << fineGridX);

  /* Right cell-left cell   pair indices: 0,1; 2,3;   4,5; 6;7
   * Front cell-back cell   pair indices: 0,2; 1,3;   4,6; 5;7
   * Top   cell-bottom cell pair indices: 0,4; 1,5;   2,6; 3;7
   *
   * Note that from the viewpoint of a cell, the face
   * has always the "opposite" index, i.e., we solve a Riemann
   * problem on the left face of the right cell (which
   * is the right face of the left cell).
   */
  constexpr int cellIndicesLeft[4] = {0, 2, 4, 6};
  constexpr int cellIndicesRight[4] = {1, 3, 5, 7};
  constexpr int cellIndicesFront[4] = {0, 1, 4, 5};
  constexpr int cellIndicesBack[4] = {2, 3, 6, 7};
#if DIMENSIONS == 3
  constexpr int cellIndicesBottom[4] = {0, 1, 2, 3};
  constexpr int cellIndicesTop[4] = {4, 5, 6, 7};
#endif
  for (int i = 0; i < TWO_POWER_D_DIVIDED_BY_TWO; i++) {
    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesLeft[i]],
      adjacentCellDescriptionsIndices[cellIndicesRight[i]],
      EXAHYPE_FACE_RIGHT, EXAHYPE_FACE_LEFT, 0);

    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesFront[i]],
      adjacentCellDescriptionsIndices[cellIndicesBack[i]],
      EXAHYPE_FACE_BACK, EXAHYPE_FACE_FRONT, 0);

    #if DIMENSIONS == 3
    solveRiemannProblemAtInterface(
      adjacentCellDescriptionsIndices[cellIndicesBottom[i]],
      adjacentCellDescriptionsIndices[cellIndicesTop[i]],
      EXAHYPE_FACE_TOP, EXAHYPE_FACE_BOTTOM, 0);
    #endif
  }
  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

void exahype::mappings::RiemannSolver::solveRiemannProblemAtInterface(
    const int cellDescriptionIndexOfLeftCell,
    const int cellDescriptionIndexOfRightCell,
    const int faceIndexForLeftCell,
    const int faceIndexForRightCell,
    const int normalNonZero) {
  // Only continue if this is an internal face, i.e.,
  // both cell description indices are valid
  if (
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfLeftCell) &&
    ADERDGCellDescriptionHeap::getInstance().isValidIndex(cellDescriptionIndexOfRightCell)
  ) {
    logDebug("touchVertexLastTime(...)::solveRiemannProblemAtInterface(...)",
             "Performing Riemann solve. "
                 << "faceIndexForLeftCell:" << faceIndexForLeftCell
                 << " faceIndexForRightCell:" << faceIndexForRightCell
                 << " indexOfLeftCell:"
                 << adjacentADERDGCellDescriptionsIndices[indexOfLeftCell]
                 << " indexOfRightCell:"
                 << adjacentADERDGCellDescriptionsIndices[indexOfRightCell]);

    std::vector<records::ADERDGCellDescription>& cellDescriptionsOfLeftCell =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfLeftCell);
    std::vector<records::ADERDGCellDescription>& cellDescriptionsOfRightCell =
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfRightCell);

    // @todo 08/02/16:Dominic Etienne Charrier
    // Assumes that the both cells hold the same number of cell descriptions
    assertion1WithExplanation(
        cellDescriptionsOfLeftCell.size() == cellDescriptionsOfRightCell.size(),
        cellDescriptionsOfLeftCell.size(),
        "The number of ADERDGCellDescriptions is not the same for both cells!");

    const int numberOfADERDGCellDescriptions = static_cast<int>(
        ADERDGCellDescriptionHeap::getInstance().getData(cellDescriptionIndexOfLeftCell).size());
    const peano::datatraversal::autotuning::MethodTrace methodTrace =
        peano::datatraversal::autotuning::UserDefined4;
    const int grainSize =
        peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
            numberOfADERDGCellDescriptions, methodTrace);

    // clang-format off
    pfor(i, 0, numberOfADERDGCellDescriptions,grainSize)
      if (cellDescriptionsOfLeftCell[i].getType() ==
                              exahype::records::ADERDGCellDescription::Cell ||
                          cellDescriptionsOfRightCell[i].getType() ==
                              exahype::records::ADERDGCellDescription::Cell) {
        exahype::solvers::Solver* solver =
            exahype::solvers::RegisteredSolvers[cellDescriptionsOfLeftCell[i]
                                                    .getSolverNumber()];

        // Lock the critical multithreading area.
        bool riemannSolveNotPerformed = false;
        assertionEquals4(
            cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell),
            cellDescriptionsOfRightCell[i].getRiemannSolvePerformed(faceIndexForRightCell),
            faceIndexForLeftCell, faceIndexForRightCell,
            cellDescriptionsOfLeftCell[i].toString(), cellDescriptionsOfRightCell[i].toString());

        riemannSolveNotPerformed =
            !cellDescriptionsOfLeftCell[i].getRiemannSolvePerformed(faceIndexForLeftCell);
        if (riemannSolveNotPerformed) {
          cellDescriptionsOfLeftCell[i].setRiemannSolvePerformed( faceIndexForLeftCell,  true);
          cellDescriptionsOfRightCell[i].setRiemannSolvePerformed(faceIndexForRightCell, true);
        }

        if (riemannSolveNotPerformed) {
          const int numberOfFaceDof = solver->getUnknownsPerFace();

          double* QL =
              DataHeap::getInstance()
                  .getData(cellDescriptionsOfLeftCell[i].getExtrapolatedPredictor())
                  .data() +
              (faceIndexForLeftCell * numberOfFaceDof);
          double* QR =
              DataHeap::getInstance()
                  .getData(cellDescriptionsOfRightCell[i].getExtrapolatedPredictor())
                  .data() +
              (faceIndexForRightCell * numberOfFaceDof);
          double* FL = DataHeap::getInstance()
                           .getData(cellDescriptionsOfLeftCell[i].getFluctuation())
                           .data() +
                       (faceIndexForLeftCell * numberOfFaceDof);
          double* FR = DataHeap::getInstance()
                           .getData(cellDescriptionsOfRightCell[i].getFluctuation())
                           .data() +
                       (faceIndexForRightCell * numberOfFaceDof);

          solver->synchroniseTimeStepping(cellDescriptionsOfLeftCell[i]);
          solver->synchroniseTimeStepping(cellDescriptionsOfRightCell[i]);

          logDebug("touchVertexLastTime(...)::debug::before::QL[0]*", QL[0]);
          logDebug("touchVertexLastTime(...)::debug::before::QR[0]*", QR[0]);
          logDebug("touchVertexLastTime(...)::debug::before::FL[0]", FL[0]);
          logDebug("touchVertexLastTime(...)::debug::before::FR[0]", FR[0]);

          solver->riemannSolver(
              FL, FR, QL, QR,
              std::min(cellDescriptionsOfLeftCell[i].getCorrectorTimeStepSize(),
                       cellDescriptionsOfRightCell[i].getCorrectorTimeStepSize()),
              normalNonZero);

          logDebug("touchVertexLastTime(...)::debug::after::QL[0]*", QL[0]);
          logDebug("touchVertexLastTime(...)::debug::after::QR[0]*", QR[0]);
          logDebug("touchVertexLastTime(...)::debug::after::FL[0]", FL[0]);
          logDebug("touchVertexLastTime(...)::debug::after::FR[0]", FR[0]);
        }
      }
    endpfor

    peano::datatraversal::autotuning::Oracle::getInstance().parallelSectionHasTerminated(methodTrace);
    // clang-format on
  }
}

void exahype::mappings::RiemannSolver::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::RiemannSolver::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  _localState = solverState;

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::RiemannSolver::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::RiemannSolver::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::RiemannSolver::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
