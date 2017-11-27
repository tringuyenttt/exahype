#include "exahype/adapters/BroadcastGlobalDataAndDropNeighbourMessages.h"


peano::CommunicationSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::communicationSpecification() const {
  return peano::CommunicationSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.communicationSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::touchVertexLastTimeSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.touchVertexLastTimeSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::touchVertexFirstTimeSpecification(int level) const { 
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.touchVertexFirstTimeSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::enterCellSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.enterCellSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::leaveCellSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.leaveCellSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::ascendSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.ascendSpecification(level)

  ;
}


peano::MappingSpecification   exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::descendSpecification(int level) const {
  return peano::MappingSpecification::getMinimalSpecification()
    &  _map2BroadcastGlobalDataAndDropNeighbourMessages.descendSpecification(level)

  ;
}


exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::BroadcastGlobalDataAndDropNeighbourMessages() {
}


exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::~BroadcastGlobalDataAndDropNeighbourMessages() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::BroadcastGlobalDataAndDropNeighbourMessages(const BroadcastGlobalDataAndDropNeighbourMessages&  masterThread):
  _map2BroadcastGlobalDataAndDropNeighbourMessages(masterThread._map2BroadcastGlobalDataAndDropNeighbourMessages) 

{
}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithWorkerThread(const BroadcastGlobalDataAndDropNeighbourMessages& workerThread) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithWorkerThread(workerThread._map2BroadcastGlobalDataAndDropNeighbourMessages);

}
#endif


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


#ifdef Parallel
void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.prepareSendToNeighbour( vertex, toRank, x, h, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}


bool exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2BroadcastGlobalDataAndDropNeighbourMessages.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithMaster(
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
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::receiveDataFromMaster(
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
   _map2BroadcastGlobalDataAndDropNeighbourMessages.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2BroadcastGlobalDataAndDropNeighbourMessages.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::beginIteration(
  exahype::State&  solverState
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.beginIteration( solverState );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::endIteration(
  exahype::State&  solverState
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.endIteration( solverState );

}




void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void exahype::adapters::BroadcastGlobalDataAndDropNeighbourMessages::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  _map2BroadcastGlobalDataAndDropNeighbourMessages.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
