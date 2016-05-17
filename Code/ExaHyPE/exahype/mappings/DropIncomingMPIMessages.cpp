#include "exahype/mappings/DropIncomingMPIMessages.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/Solver.h"


peano::CommunicationSpecification   exahype::mappings::DropIncomingMPIMessages::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true
  );
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid);
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces);
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


peano::MappingSpecification   exahype::mappings::DropIncomingMPIMessages::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidCoarseGridRaces);
}


tarch::logging::Log                exahype::mappings::DropIncomingMPIMessages::_log( "exahype::mappings::DropIncomingMPIMessages" ); 


exahype::mappings::DropIncomingMPIMessages::DropIncomingMPIMessages() {
}


exahype::mappings::DropIncomingMPIMessages::~DropIncomingMPIMessages() {
}


#if defined(SharedMemoryParallelisation)
exahype::mappings::DropIncomingMPIMessages::DropIncomingMPIMessages(const DropIncomingMPIMessages&  masterThread) {
}


void exahype::mappings::DropIncomingMPIMessages::mergeWithWorkerThread(const DropIncomingMPIMessages& workerThread) {
}
#endif


void exahype::mappings::DropIncomingMPIMessages::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}

#ifdef Parallel
void exahype::mappings::DropIncomingMPIMessages::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  logTraceInWith6Arguments( "mergeWithNeighbour(...)", vertex, neighbour, fromRank, fineGridX, fineGridH, level );

  #if !defined(PeriodicBC)
  if (vertex.isBoundary()) return;
  #endif

  tarch::la::Vector<TWO_POWER_D, int>& adjacentADERDGCellDescriptionsIndices =
      vertex.getADERDGCellDescriptionsIndex();

  dfor2(dest)
  dfor2(src)
    if (
      vertex.getAdjacentRanks()(destScalar)==tarch::parallel::Node::getInstance().getRank()
      &&
      vertex.getAdjacentRanks()(srcScalar)==fromRank
      &&
      tarch::la::countEqualEntries(dest,src)==1    // we are solely exchanging faces
    ) {
      const int destCellDescriptionIndex = adjacentADERDGCellDescriptionsIndices(destScalar);
      assertion5(
        ADERDGCellDescriptionHeap::getInstance().isValidIndex(destCellDescriptionIndex),
        src, dest,
        multiscalelinkedcell::indicesToString( adjacentADERDGCellDescriptionsIndices ),
        vertex.toString(),
        tarch::parallel::Node::getInstance().getRank()
      );
      std::vector<records::ADERDGCellDescription>& cellDescriptions = ADERDGCellDescriptionHeap::getInstance().getData(destCellDescriptionIndex);

      for (int currentSolver=0; currentSolver<static_cast<int>(cellDescriptions.size()); currentSolver++) {
        if (cellDescriptions[currentSolver].getType()==exahype::records::ADERDGCellDescription::Cell) {
          exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[cellDescriptions[currentSolver].getSolverNumber()];

          const int numberOfFaceDof       = solver->getUnknownsPerFace();
          const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src,dest);
          assertion(normalOfExchangedFace>=0 && normalOfExchangedFace<DIMENSIONS);
          const int offsetInFaceArray     = 2*normalOfExchangedFace + src(normalOfExchangedFace)<dest(normalOfExchangedFace) ? 0 : 1;

          assertion( DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getExtrapolatedPredictor()) );
          assertion( DataHeap::getInstance().isValidIndex(cellDescriptions[currentSolver].getFluctuation()) );

          const double* lQhbnd = DataHeap::getInstance().getData(cellDescriptions[currentSolver].getExtrapolatedPredictor()).data() + (offsetInFaceArray * numberOfFaceDof);
          const double* lFhbnd = DataHeap::getInstance().getData(cellDescriptions[currentSolver].getFluctuation()).data()           + (offsetInFaceArray * numberOfFaceDof);

          if ( adjacentADERDGCellDescriptionsIndices(destScalar)==multiscalelinkedcell::HangingVertexBookkeeper::DomainBoundaryAdjacencyIndex ) {
            #if defined(PeriodicBC)
            assertionMsg( false, "Vasco, we have to implement this" );
            #else
            assertionMsg( false, "should never been entered");
            #endif
          }
          else {
            logDebug( "mergeWithNeighbour(...)", "receive two arrays from rank " << fromRank << " for vertex " << vertex.toString() << ", src type=" << multiscalelinkedcell::indexToString(adjacentADERDGCellDescriptionsIndices(srcScalar)) );

            DataHeap::HeapEntries lQhbnd = DataHeap::getInstance().receiveData( fromRank, fineGridX, level, peano::heap::MessageType::NeighbourCommunication );
            DataHeap::HeapEntries lFhbnd = DataHeap::getInstance().receiveData( fromRank, fineGridX, level, peano::heap::MessageType::NeighbourCommunication );
          }
        }
        else {
          assertionMsg( false, "Dominic, please implement" );
        }
      }
    }
  enddforx
  enddforx

  logTraceOut( "mergeWithNeighbour(...)" );
}

void exahype::mappings::DropIncomingMPIMessages::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}

void exahype::mappings::DropIncomingMPIMessages::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}

void exahype::mappings::DropIncomingMPIMessages::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}

void exahype::mappings::DropIncomingMPIMessages::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}

void exahype::mappings::DropIncomingMPIMessages::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                       level
) {
}

bool exahype::mappings::DropIncomingMPIMessages::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}

void exahype::mappings::DropIncomingMPIMessages::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::mergeWithMaster(
  const exahype::Cell&           workerGridCell,
  exahype::Vertex * const        workerGridVertices,
 const peano::grid::VertexEnumerator& workerEnumerator,
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
  const exahype::State&          workerState,
  exahype::State&                masterState
) {
}


void exahype::mappings::DropIncomingMPIMessages::receiveDataFromMaster(
      exahype::Cell&                        receivedCell, 
      exahype::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      exahype::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      exahype::Cell&                        receivedCoarseGridCell,
      exahype::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      exahype::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
}


void exahype::mappings::DropIncomingMPIMessages::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
}
#endif

void exahype::mappings::DropIncomingMPIMessages::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::mappings::DropIncomingMPIMessages::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::beginIteration(
  exahype::State&  solverState
) {
}


void exahype::mappings::DropIncomingMPIMessages::endIteration(
  exahype::State&  solverState
) {
}



void exahype::mappings::DropIncomingMPIMessages::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
}


void exahype::mappings::DropIncomingMPIMessages::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
}