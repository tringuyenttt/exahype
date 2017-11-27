#include "exahype/repositories/RepositoryArrayStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#ifdef Parallel
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#endif

#include "peano/datatraversal/autotuning/Oracle.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "peano/performanceanalysis/ScorePMacros.h"

#if !defined(CompilerICC)
#include "peano/grid/Grid.cpph"
#endif


tarch::logging::Log exahype::repositories::RepositoryArrayStack::_log( "exahype::repositories::RepositoryArrayStack" );


exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack, maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithMeshRefinementAndPlotGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFinaliseMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFinaliseMeshRefinementOrLocalRollback(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFusedTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithBroadcastGlobalDataAndDropNeighbourMessages(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionOrLocalRecomputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalRollback(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithBroadcastGlobalDataAndMergeNeighbourMessages(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(...)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(...)" );
}



exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack,maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithMeshRefinementAndPlotGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFinaliseMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFinaliseMeshRefinementOrLocalRollback(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithFusedTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithBroadcastGlobalDataAndDropNeighbourMessages(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionOrLocalRecomputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGlobalRollback(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithBroadcastGlobalDataAndMergeNeighbourMessages(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(Geometry&)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(Geometry&)" );
}
    
   
exahype::repositories::RepositoryArrayStack::~RepositoryArrayStack() {
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );
}


void exahype::repositories::RepositoryArrayStack::restart(
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          domainLevel,
  const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
) {
  logTraceInWith4Arguments( "restart(...)", domainSize, domainOffset, domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel );
  #ifdef Parallel
  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());
  #endif
  
  logInfo( "restart(...)", "start node for subdomain " << domainOffset << "x" << domainSize << " on level " << domainLevel << " with master " << tarch::parallel::NodePool::getInstance().getMasterRank() );
  
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithMeshRefinement.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithMeshRefinementAndPlotGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFinaliseMeshRefinement.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFinaliseMeshRefinementOrLocalRollback.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithFusedTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionRerun.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithBroadcastGlobalDataAndDropNeighbourMessages.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithLimiterStatusSpreading.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionOrLocalRecomputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGlobalRollback.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithBroadcastGlobalDataAndMergeNeighbourMessages.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionUpdate.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPrediction.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);

 
   _solverState.restart();
 
  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositoryArrayStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  peano::parallel::SendReceiveBufferPool::getInstance().terminate();
  #endif
  
  _gridWithMeshRefinement.terminate();
  _gridWithMeshRefinementAndPlotGrid.terminate();
  _gridWithFinaliseMeshRefinement.terminate();
  _gridWithFinaliseMeshRefinementOrLocalRollback.terminate();
  _gridWithFusedTimeStep.terminate();
  _gridWithPredictionRerun.terminate();
  _gridWithBroadcastGlobalDataAndDropNeighbourMessages.terminate();
  _gridWithLimiterStatusSpreading.terminate();
  _gridWithPredictionOrLocalRecomputation.terminate();
  _gridWithGlobalRollback.terminate();
  _gridWithBroadcastGlobalDataAndMergeNeighbourMessages.terminate();
  _gridWithSolutionUpdate.terminate();
  _gridWithPrediction.terminate();

 
  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositoryArrayStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositoryArrayStack::getState() const {
  return _solverState;
}

   
void exahype::repositories::RepositoryArrayStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  SCOREP_USER_REGION( (std::string("exahype::repositories::RepositoryArrayStack::iterate() - ") + _repositoryState.toString( _repositoryState.getAction() )).c_str(), SCOREP_USER_REGION_TYPE_FUNCTION)

  tarch::timing::Watch watch( "exahype::repositories::RepositoryArrayStack", "iterate(bool)", false);
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    _repositoryState.setNumberOfIterations(numberOfIterations);
    _repositoryState.setExchangeBoundaryVertices(exchangeBoundaryVertices);
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  else {
    assertionEquals( numberOfIterations, 1 );
    numberOfIterations = _repositoryState.getNumberOfIterations();
  }

  peano::parallel::SendReceiveBufferPool::getInstance().exchangeBoundaryVertices(_repositoryState.getExchangeBoundaryVertices());

  if ( numberOfIterations > 1 && _solverState.isInvolvedInJoinOrFork() ) {
    logWarning( "iterate()", "iterate invoked for multiple traversals though load balancing still does redistribute data" );
  }
  bool switchedLoadBalancingTemporarilyOff = false;
  if ( numberOfIterations > 1 && peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated() ) {
    switchedLoadBalancingTemporarilyOff = true;
    peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(false);
  }

  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());

  peano::parallel::loadbalancing::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {

    #ifdef Parallel
    if (_repositoryState.getNumberOfIterations()==1) {
      _solverState.currentlyRunsMultipleIterations(State::BatchState::NoBatch);
    }
    else if (i==0) {
      _solverState.currentlyRunsMultipleIterations(State::BatchState::FirstIterationOfBatch);
    }
    else if (i==numberOfIterations-1) {
      _solverState.currentlyRunsMultipleIterations(State::BatchState::LastIterationOfBatch);
    }
    else {
      _solverState.currentlyRunsMultipleIterations(State::BatchState::IntermediateIterationOfBatch);
    }
    #endif
  
  
    switch ( _repositoryState.getAction()) {
      case exahype::records::RepositoryState::UseAdapterMeshRefinement: watch.startTimer(); _gridWithMeshRefinement.iterate(); watch.stopTimer(); _measureMeshRefinementCPUTime.setValue( watch.getCPUTime() ); _measureMeshRefinementCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterMeshRefinementAndPlotGrid: watch.startTimer(); _gridWithMeshRefinementAndPlotGrid.iterate(); watch.stopTimer(); _measureMeshRefinementAndPlotGridCPUTime.setValue( watch.getCPUTime() ); _measureMeshRefinementAndPlotGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinement: watch.startTimer(); _gridWithFinaliseMeshRefinement.iterate(); watch.stopTimer(); _measureFinaliseMeshRefinementCPUTime.setValue( watch.getCPUTime() ); _measureFinaliseMeshRefinementCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinementOrLocalRollback: watch.startTimer(); _gridWithFinaliseMeshRefinementOrLocalRollback.iterate(); watch.stopTimer(); _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.setValue( watch.getCPUTime() ); _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterFusedTimeStep: watch.startTimer(); _gridWithFusedTimeStep.iterate(); watch.stopTimer(); _measureFusedTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureFusedTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionRerun: watch.startTimer(); _gridWithPredictionRerun.iterate(); watch.stopTimer(); _measurePredictionRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictionRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndDropNeighbourMessages: watch.startTimer(); _gridWithBroadcastGlobalDataAndDropNeighbourMessages.iterate(); watch.stopTimer(); _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.setValue( watch.getCPUTime() ); _measureBroadcastGlobalDataAndDropNeighbourMessagesCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading: watch.startTimer(); _gridWithLimiterStatusSpreading.iterate(); watch.stopTimer(); _measureLimiterStatusSpreadingCPUTime.setValue( watch.getCPUTime() ); _measureLimiterStatusSpreadingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionOrLocalRecomputation: watch.startTimer(); _gridWithPredictionOrLocalRecomputation.iterate(); watch.stopTimer(); _measurePredictionOrLocalRecomputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionOrLocalRecomputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGlobalRollback: watch.startTimer(); _gridWithGlobalRollback.iterate(); watch.stopTimer(); _measureGlobalRollbackCPUTime.setValue( watch.getCPUTime() ); _measureGlobalRollbackCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndMergeNeighbourMessages: watch.startTimer(); _gridWithBroadcastGlobalDataAndMergeNeighbourMessages.iterate(); watch.stopTimer(); _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.setValue( watch.getCPUTime() ); _measureBroadcastGlobalDataAndMergeNeighbourMessagesCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionUpdate: watch.startTimer(); _gridWithSolutionUpdate.iterate(); watch.stopTimer(); _measureSolutionUpdateCPUTime.setValue( watch.getCPUTime() ); _measureSolutionUpdateCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPrediction: watch.startTimer(); _gridWithPrediction.iterate(); watch.stopTimer(); _measurePredictionCPUTime.setValue( watch.getCPUTime() ); _measurePredictionCalendarTime.setValue( watch.getCalendarTime() ); break;

      case exahype::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case exahype::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
    #ifdef Parallel
    if ( switchedLoadBalancingTemporarilyOff && i==numberOfIterations-1) {
      peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(true);
    }
    #endif
  }
    
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  }
  #endif
}

 void exahype::repositories::RepositoryArrayStack::switchToMeshRefinement() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterMeshRefinement); }
 void exahype::repositories::RepositoryArrayStack::switchToMeshRefinementAndPlotGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterMeshRefinementAndPlotGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToFinaliseMeshRefinement() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinement); }
 void exahype::repositories::RepositoryArrayStack::switchToFinaliseMeshRefinementOrLocalRollback() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinementOrLocalRollback); }
 void exahype::repositories::RepositoryArrayStack::switchToFusedTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterFusedTimeStep); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionRerun); }
 void exahype::repositories::RepositoryArrayStack::switchToBroadcastGlobalDataAndDropNeighbourMessages() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndDropNeighbourMessages); }
 void exahype::repositories::RepositoryArrayStack::switchToLimiterStatusSpreading() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionOrLocalRecomputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionOrLocalRecomputation); }
 void exahype::repositories::RepositoryArrayStack::switchToGlobalRollback() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGlobalRollback); }
 void exahype::repositories::RepositoryArrayStack::switchToBroadcastGlobalDataAndMergeNeighbourMessages() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndMergeNeighbourMessages); }
 void exahype::repositories::RepositoryArrayStack::switchToSolutionUpdate() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionUpdate); }
 void exahype::repositories::RepositoryArrayStack::switchToPrediction() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPrediction); }



 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterMeshRefinement() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterMeshRefinement; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterMeshRefinementAndPlotGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterMeshRefinementAndPlotGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterFinaliseMeshRefinement() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinement; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterFinaliseMeshRefinementOrLocalRollback() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFinaliseMeshRefinementOrLocalRollback; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterFusedTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterFusedTimeStep; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionRerun; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterBroadcastGlobalDataAndDropNeighbourMessages() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndDropNeighbourMessages; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterLimiterStatusSpreading() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionOrLocalRecomputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionOrLocalRecomputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterGlobalRollback() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGlobalRollback; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterBroadcastGlobalDataAndMergeNeighbourMessages() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterBroadcastGlobalDataAndMergeNeighbourMessages; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterSolutionUpdate() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionUpdate; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPrediction() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPrediction; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositoryArrayStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositoryArrayStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositoryArrayStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


void exahype::repositories::RepositoryArrayStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


#ifdef Parallel
void exahype::repositories::RepositoryArrayStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositoryArrayStack::ContinueCommand exahype::repositories::RepositoryArrayStack::continueToIterate() {
  logTraceIn( "continueToIterate()" );

  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());

  ContinueCommand result;
  if ( _solverState.hasJoinedWithMaster() ) {
    result = Terminate;
  }
  else {
    int masterNode = tarch::parallel::Node::getInstance().getGlobalMasterRank();
    assertion( masterNode != -1 );

    _repositoryState.receive( masterNode, peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag(), true, ReceiveIterationControlMessagesBlocking );

    result = Continue;
    if (_repositoryState.getAction()==exahype::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==exahype::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void exahype::repositories::RepositoryArrayStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureMeshRefinementCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| MeshRefinement \t |  " << _measureMeshRefinementCPUTime.getNumberOfMeasurements() << " \t |  " << _measureMeshRefinementCPUTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCPUTime.getValue()  << " \t |  " << _measureMeshRefinementCalendarTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCalendarTime.getValue() << " \t |  " << _measureMeshRefinementCPUTime.toString() << " \t |  " << _measureMeshRefinementCalendarTime.toString() );
   if (logAllAdapters || _measureMeshRefinementAndPlotGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| MeshRefinementAndPlotGrid \t |  " << _measureMeshRefinementAndPlotGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureMeshRefinementAndPlotGridCPUTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementAndPlotGridCPUTime.getValue()  << " \t |  " << _measureMeshRefinementAndPlotGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementAndPlotGridCalendarTime.getValue() << " \t |  " << _measureMeshRefinementAndPlotGridCPUTime.toString() << " \t |  " << _measureMeshRefinementAndPlotGridCalendarTime.toString() );
   if (logAllAdapters || _measureFinaliseMeshRefinementCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| FinaliseMeshRefinement \t |  " << _measureFinaliseMeshRefinementCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFinaliseMeshRefinementCPUTime.getAccumulatedValue() << " \t |  " << _measureFinaliseMeshRefinementCPUTime.getValue()  << " \t |  " << _measureFinaliseMeshRefinementCalendarTime.getAccumulatedValue() << " \t |  " << _measureFinaliseMeshRefinementCalendarTime.getValue() << " \t |  " << _measureFinaliseMeshRefinementCPUTime.toString() << " \t |  " << _measureFinaliseMeshRefinementCalendarTime.toString() );
   if (logAllAdapters || _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| FinaliseMeshRefinementOrLocalRollback \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.getAccumulatedValue() << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.getValue()  << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime.getAccumulatedValue() << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime.getValue() << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.toString() << " \t |  " << _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime.toString() );
   if (logAllAdapters || _measureFusedTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| FusedTimeStep \t |  " << _measureFusedTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureFusedTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureFusedTimeStepCPUTime.getValue()  << " \t |  " << _measureFusedTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureFusedTimeStepCalendarTime.getValue() << " \t |  " << _measureFusedTimeStepCPUTime.toString() << " \t |  " << _measureFusedTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionRerun \t |  " << _measurePredictionRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCPUTime.getValue()  << " \t |  " << _measurePredictionRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCalendarTime.getValue() << " \t |  " << _measurePredictionRerunCPUTime.toString() << " \t |  " << _measurePredictionRerunCalendarTime.toString() );
   if (logAllAdapters || _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| BroadcastGlobalDataAndDropNeighbourMessages \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.getNumberOfMeasurements() << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.getAccumulatedValue() << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.getValue()  << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCalendarTime.getAccumulatedValue() << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCalendarTime.getValue() << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.toString() << " \t |  " << _measureBroadcastGlobalDataAndDropNeighbourMessagesCalendarTime.toString() );
   if (logAllAdapters || _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| LimiterStatusSpreading \t |  " << _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getValue()  << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.toString() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionOrLocalRecomputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionOrLocalRecomputation \t |  " << _measurePredictionOrLocalRecomputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionOrLocalRecomputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionOrLocalRecomputationCPUTime.getValue()  << " \t |  " << _measurePredictionOrLocalRecomputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionOrLocalRecomputationCalendarTime.getValue() << " \t |  " << _measurePredictionOrLocalRecomputationCPUTime.toString() << " \t |  " << _measurePredictionOrLocalRecomputationCalendarTime.toString() );
   if (logAllAdapters || _measureGlobalRollbackCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GlobalRollback \t |  " << _measureGlobalRollbackCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGlobalRollbackCPUTime.getAccumulatedValue() << " \t |  " << _measureGlobalRollbackCPUTime.getValue()  << " \t |  " << _measureGlobalRollbackCalendarTime.getAccumulatedValue() << " \t |  " << _measureGlobalRollbackCalendarTime.getValue() << " \t |  " << _measureGlobalRollbackCPUTime.toString() << " \t |  " << _measureGlobalRollbackCalendarTime.toString() );
   if (logAllAdapters || _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| BroadcastGlobalDataAndMergeNeighbourMessages \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.getNumberOfMeasurements() << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.getAccumulatedValue() << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.getValue()  << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCalendarTime.getAccumulatedValue() << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCalendarTime.getValue() << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.toString() << " \t |  " << _measureBroadcastGlobalDataAndMergeNeighbourMessagesCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionUpdateCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionUpdate \t |  " << _measureSolutionUpdateCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionUpdateCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCPUTime.getValue()  << " \t |  " << _measureSolutionUpdateCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCalendarTime.getValue() << " \t |  " << _measureSolutionUpdateCPUTime.toString() << " \t |  " << _measureSolutionUpdateCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Prediction \t |  " << _measurePredictionCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionCPUTime.getValue()  << " \t |  " << _measurePredictionCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionCalendarTime.getValue() << " \t |  " << _measurePredictionCPUTime.toString() << " \t |  " << _measurePredictionCalendarTime.toString() );

}


void exahype::repositories::RepositoryArrayStack::clearIterationStatistics() {
   _measureMeshRefinementCPUTime.erase();
   _measureMeshRefinementAndPlotGridCPUTime.erase();
   _measureFinaliseMeshRefinementCPUTime.erase();
   _measureFinaliseMeshRefinementOrLocalRollbackCPUTime.erase();
   _measureFusedTimeStepCPUTime.erase();
   _measurePredictionRerunCPUTime.erase();
   _measureBroadcastGlobalDataAndDropNeighbourMessagesCPUTime.erase();
   _measureLimiterStatusSpreadingCPUTime.erase();
   _measurePredictionOrLocalRecomputationCPUTime.erase();
   _measureGlobalRollbackCPUTime.erase();
   _measureBroadcastGlobalDataAndMergeNeighbourMessagesCPUTime.erase();
   _measureSolutionUpdateCPUTime.erase();
   _measurePredictionCPUTime.erase();

   _measureMeshRefinementCalendarTime.erase();
   _measureMeshRefinementAndPlotGridCalendarTime.erase();
   _measureFinaliseMeshRefinementCalendarTime.erase();
   _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime.erase();
   _measureFusedTimeStepCalendarTime.erase();
   _measurePredictionRerunCalendarTime.erase();
   _measureBroadcastGlobalDataAndDropNeighbourMessagesCalendarTime.erase();
   _measureLimiterStatusSpreadingCalendarTime.erase();
   _measurePredictionOrLocalRecomputationCalendarTime.erase();
   _measureGlobalRollbackCalendarTime.erase();
   _measureBroadcastGlobalDataAndMergeNeighbourMessagesCalendarTime.erase();
   _measureSolutionUpdateCalendarTime.erase();
   _measurePredictionCalendarTime.erase();

}
