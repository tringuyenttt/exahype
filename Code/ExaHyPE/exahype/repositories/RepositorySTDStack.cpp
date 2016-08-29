#include "exahype/repositories/RepositorySTDStack.h"

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

#if !defined(CompilerICC)
#include "peano/grid/Grid.cpph"
#endif


tarch::logging::Log exahype::repositories::RepositorySTDStack::_log( "exahype::repositories::RepositorySTDStack" );


exahype::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(...)" );
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(...)" );
}



exahype::repositories::RepositorySTDStack::RepositorySTDStack(
  peano::geometry::Geometry&                   geometry
):
  _geometry(geometry),
  _cellStack(),
  _vertexStack(),
  _solverState(),
  _gridWithAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorAndGlobalTimeStepComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStepAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictorRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithRiemannSolver(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictor(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrector(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithCorrectorAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositorySTDStack(Geometry&)" );

  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositorySTDStack(Geometry&)" );
}
    
   
exahype::repositories::RepositorySTDStack::~RepositorySTDStack() {
  assertionMsg( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate, "terminate() must be called before destroying repository." );
}


void exahype::repositories::RepositorySTDStack::restart(
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

  _gridWithAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorAndGlobalTimeStepComputation.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGridErasing.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStepAndPlot.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictorRerun.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithRiemannSolver.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictor.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrector.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithCorrectorAndPlot.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlot.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);


  _solverState.restart();

  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositorySTDStack::terminate() {
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

  _gridWithAugmentedAMRGrid.terminate();
  _gridWithPlotAugmentedAMRGrid.terminate();
  _gridWithInitialConditionAndGlobalTimeStepComputation.terminate();
  _gridWithPredictorAndPlotAndGlobalTimeStepComputation.terminate();
  _gridWithPredictorAndGlobalTimeStepComputation.terminate();
  _gridWithGridErasing.terminate();
  _gridWithADERDGTimeStep.terminate();
  _gridWithADERDGTimeStepAndPlot.terminate();
  _gridWithPredictorRerun.terminate();
  _gridWithRiemannSolver.terminate();
  _gridWithPredictor.terminate();
  _gridWithCorrector.terminate();
  _gridWithCorrectorAndPlot.terminate();
  _gridWithPlot.terminate();


  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositorySTDStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositorySTDStack::getState() const {
  return _solverState;
}


void exahype::repositories::RepositorySTDStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "exahype::repositories::RepositorySTDStack", "iterate(bool)", false);
  
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
  
  _solverState.currentlyRunsMultipleIterations(_repositoryState.getNumberOfIterations()>1);
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {
    switch ( _repositoryState.getAction()) {
      case exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid: watch.startTimer(); _gridWithAugmentedAMRGrid.iterate(); watch.stopTimer(); _measureAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measureAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid: watch.startTimer(); _gridWithPlotAugmentedAMRGrid.iterate(); watch.stopTimer(); _measurePlotAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measurePlotAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation: watch.startTimer(); _gridWithInitialConditionAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPredictorAndPlotAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation: watch.startTimer(); _gridWithPredictorAndGlobalTimeStepComputation.iterate(); watch.stopTimer(); _measurePredictorAndGlobalTimeStepComputationCPUTime.setValue( watch.getCPUTime() ); _measurePredictorAndGlobalTimeStepComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGridErasing: watch.startTimer(); _gridWithGridErasing.iterate(); watch.stopTimer(); _measureGridErasingCPUTime.setValue( watch.getCPUTime() ); _measureGridErasingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStep: watch.startTimer(); _gridWithADERDGTimeStep.iterate(); watch.stopTimer(); _measureADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot: watch.startTimer(); _gridWithADERDGTimeStepAndPlot.iterate(); watch.stopTimer(); _measureADERDGTimeStepAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictorRerun: watch.startTimer(); _gridWithPredictorRerun.iterate(); watch.stopTimer(); _measurePredictorRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictorRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterRiemannSolver: watch.startTimer(); _gridWithRiemannSolver.iterate(); watch.stopTimer(); _measureRiemannSolverCPUTime.setValue( watch.getCPUTime() ); _measureRiemannSolverCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictor: watch.startTimer(); _gridWithPredictor.iterate(); watch.stopTimer(); _measurePredictorCPUTime.setValue( watch.getCPUTime() ); _measurePredictorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrector: watch.startTimer(); _gridWithCorrector.iterate(); watch.stopTimer(); _measureCorrectorCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterCorrectorAndPlot: watch.startTimer(); _gridWithCorrectorAndPlot.iterate(); watch.stopTimer(); _measureCorrectorAndPlotCPUTime.setValue( watch.getCPUTime() ); _measureCorrectorAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlot: watch.startTimer(); _gridWithPlot.iterate(); watch.stopTimer(); _measurePlotCPUTime.setValue( watch.getCPUTime() ); _measurePlotCalendarTime.setValue( watch.getCalendarTime() ); break;

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

 void exahype::repositories::RepositorySTDStack::switchToAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid); }
 void exahype::repositories::RepositorySTDStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositorySTDStack::switchToInitialConditionAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorAndPlotAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorAndGlobalTimeStepComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation); }
 void exahype::repositories::RepositorySTDStack::switchToGridErasing() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGridErasing); }
 void exahype::repositories::RepositorySTDStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositorySTDStack::switchToADERDGTimeStepAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot); }
 void exahype::repositories::RepositorySTDStack::switchToPredictorRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictorRerun); }
 void exahype::repositories::RepositorySTDStack::switchToRiemannSolver() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterRiemannSolver); }
 void exahype::repositories::RepositorySTDStack::switchToPredictor() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictor); }
 void exahype::repositories::RepositorySTDStack::switchToCorrector() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrector); }
 void exahype::repositories::RepositorySTDStack::switchToCorrectorAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrectorAndPlot); }
 void exahype::repositories::RepositorySTDStack::switchToPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlot); }



 bool exahype::repositories::RepositorySTDStack::isActiveAdapterAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterAugmentedAMRGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterInitialConditionAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorAndPlotAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndPlotAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorAndGlobalTimeStepComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorAndGlobalTimeStepComputation; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterGridErasing() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridErasing; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterADERDGTimeStepAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStepAndPlot; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictorRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictorRerun; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterRiemannSolver() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterRiemannSolver; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPredictor() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictor; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrector() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrector; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrectorAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrectorAndPlot; }
 bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlot; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositorySTDStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositorySTDStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositorySTDStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


void exahype::repositories::RepositorySTDStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


#ifdef Parallel
void exahype::repositories::RepositorySTDStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositorySTDStack::ContinueCommand exahype::repositories::RepositorySTDStack::continueToIterate() {
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


void exahype::repositories::RepositorySTDStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| AugmentedAMRGrid \t |  " << _measureAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measureAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measureAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measureAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measureAugmentedAMRGridCPUTime.toString() << " \t |  " << _measureAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| InitialConditionAndGlobalTimeStepComputation \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndPlotAndGlobalTimeStepComputation \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorAndGlobalTimeStepComputation \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.getValue()  << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.getValue() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCPUTime.toString() << " \t |  " << _measurePredictorAndGlobalTimeStepComputationCalendarTime.toString() );
   if (logAllAdapters || _measureGridErasingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GridErasing \t |  " << _measureGridErasingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGridErasingCPUTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCPUTime.getValue()  << " \t |  " << _measureGridErasingCalendarTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCalendarTime.getValue() << " \t |  " << _measureGridErasingCPUTime.toString() << " \t |  " << _measureGridErasingCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStepAndPlot \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepAndPlotCPUTime.toString() << " \t |  " << _measureADERDGTimeStepAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictorRerun \t |  " << _measurePredictorRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCPUTime.getValue()  << " \t |  " << _measurePredictorRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorRerunCalendarTime.getValue() << " \t |  " << _measurePredictorRerunCPUTime.toString() << " \t |  " << _measurePredictorRerunCalendarTime.toString() );
   if (logAllAdapters || _measureRiemannSolverCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| RiemannSolver \t |  " << _measureRiemannSolverCPUTime.getNumberOfMeasurements() << " \t |  " << _measureRiemannSolverCPUTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCPUTime.getValue()  << " \t |  " << _measureRiemannSolverCalendarTime.getAccumulatedValue() << " \t |  " << _measureRiemannSolverCalendarTime.getValue() << " \t |  " << _measureRiemannSolverCPUTime.toString() << " \t |  " << _measureRiemannSolverCalendarTime.toString() );
   if (logAllAdapters || _measurePredictorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Predictor \t |  " << _measurePredictorCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictorCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictorCPUTime.getValue()  << " \t |  " << _measurePredictorCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictorCalendarTime.getValue() << " \t |  " << _measurePredictorCPUTime.toString() << " \t |  " << _measurePredictorCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectorCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Corrector \t |  " << _measureCorrectorCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCPUTime.getValue()  << " \t |  " << _measureCorrectorCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorCalendarTime.getValue() << " \t |  " << _measureCorrectorCPUTime.toString() << " \t |  " << _measureCorrectorCalendarTime.toString() );
   if (logAllAdapters || _measureCorrectorAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| CorrectorAndPlot \t |  " << _measureCorrectorAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectorAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPlotCPUTime.getValue()  << " \t |  " << _measureCorrectorAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectorAndPlotCalendarTime.getValue() << " \t |  " << _measureCorrectorAndPlotCPUTime.toString() << " \t |  " << _measureCorrectorAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Plot \t |  " << _measurePlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotCPUTime.getValue()  << " \t |  " << _measurePlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotCalendarTime.getValue() << " \t |  " << _measurePlotCPUTime.toString() << " \t |  " << _measurePlotCalendarTime.toString() );

}


void exahype::repositories::RepositorySTDStack::clearIterationStatistics() {
   _measureAugmentedAMRGridCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureInitialConditionAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCPUTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCPUTime.erase();
   _measureGridErasingCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measureADERDGTimeStepAndPlotCPUTime.erase();
   _measurePredictorRerunCPUTime.erase();
   _measureRiemannSolverCPUTime.erase();
   _measurePredictorCPUTime.erase();
   _measureCorrectorCPUTime.erase();
   _measureCorrectorAndPlotCPUTime.erase();
   _measurePlotCPUTime.erase();

   _measureAugmentedAMRGridCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureInitialConditionAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndPlotAndGlobalTimeStepComputationCalendarTime.erase();
   _measurePredictorAndGlobalTimeStepComputationCalendarTime.erase();
   _measureGridErasingCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measureADERDGTimeStepAndPlotCalendarTime.erase();
   _measurePredictorRerunCalendarTime.erase();
   _measureRiemannSolverCalendarTime.erase();
   _measurePredictorCalendarTime.erase();
   _measureCorrectorCalendarTime.erase();
   _measureCorrectorAndPlotCalendarTime.erase();
   _measurePlotCalendarTime.erase();

}
