// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STACK_H_ 


#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/VertexArrayStack.h"


 #include "exahype/adapters/MeshRefinement.h" 
 #include "exahype/adapters/GridErasing.h" 
 #include "exahype/adapters/FusedTimeStep.h" 
 #include "exahype/adapters/PlotAndFusedTimeStep.h" 
 #include "exahype/adapters/LimiterStatusSpreading.h" 
 #include "exahype/adapters/Reinitialisation.h" 
 #include "exahype/adapters/LocalRecomputationAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/GlobalRollback.h" 
 #include "exahype/adapters/NeighbourDataMerging.h" 
 #include "exahype/adapters/SolutionUpdateAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/TimeStepSizeComputation.h" 
 #include "exahype/adapters/Prediction.h" 
 #include "exahype/adapters/PredictionAndPlot.h" 
 #include "exahype/adapters/FinaliseMeshRefinementAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/MergeTimeStepData.h" 
 #include "exahype/adapters/MergeTimeStepDataDropFaceData.h" 
 #include "exahype/adapters/FinaliseMeshRefinementAndReinitialisation.h" 



namespace exahype {
      namespace repositories {
        class RepositoryArrayStack;  
      }
}


class exahype::repositories::RepositoryArrayStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellArrayStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexArrayStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MeshRefinement> _gridWithMeshRefinement;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GridErasing> _gridWithGridErasing;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FusedTimeStep> _gridWithFusedTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAndFusedTimeStep> _gridWithPlotAndFusedTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusSpreading> _gridWithLimiterStatusSpreading;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Reinitialisation> _gridWithReinitialisation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LocalRecomputationAndTimeStepSizeComputation> _gridWithLocalRecomputationAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalRollback> _gridWithGlobalRollback;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::NeighbourDataMerging> _gridWithNeighbourDataMerging;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionUpdateAndTimeStepSizeComputation> _gridWithSolutionUpdateAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::TimeStepSizeComputation> _gridWithTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Prediction> _gridWithPrediction;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndPlot> _gridWithPredictionAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinementAndTimeStepSizeComputation> _gridWithFinaliseMeshRefinementAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MergeTimeStepData> _gridWithMergeTimeStepData;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MergeTimeStepDataDropFaceData> _gridWithMergeTimeStepDataDropFaceData;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinementAndReinitialisation> _gridWithFinaliseMeshRefinementAndReinitialisation;

  
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureMeshRefinementCPUTime;
    tarch::timing::Measurement _measureGridErasingCPUTime;
    tarch::timing::Measurement _measureFusedTimeStepCPUTime;
    tarch::timing::Measurement _measurePlotAndFusedTimeStepCPUTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCPUTime;
    tarch::timing::Measurement _measureReinitialisationCPUTime;
    tarch::timing::Measurement _measureLocalRecomputationAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measureGlobalRollbackCPUTime;
    tarch::timing::Measurement _measureNeighbourDataMergingCPUTime;
    tarch::timing::Measurement _measureSolutionUpdateAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measureTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measurePredictionCPUTime;
    tarch::timing::Measurement _measurePredictionAndPlotCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measureMergeTimeStepDataCPUTime;
    tarch::timing::Measurement _measureMergeTimeStepDataDropFaceDataCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementAndReinitialisationCPUTime;

    tarch::timing::Measurement _measureMeshRefinementCalendarTime;
    tarch::timing::Measurement _measureGridErasingCalendarTime;
    tarch::timing::Measurement _measureFusedTimeStepCalendarTime;
    tarch::timing::Measurement _measurePlotAndFusedTimeStepCalendarTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCalendarTime;
    tarch::timing::Measurement _measureReinitialisationCalendarTime;
    tarch::timing::Measurement _measureLocalRecomputationAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measureGlobalRollbackCalendarTime;
    tarch::timing::Measurement _measureNeighbourDataMergingCalendarTime;
    tarch::timing::Measurement _measureSolutionUpdateAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measureTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measurePredictionCalendarTime;
    tarch::timing::Measurement _measurePredictionAndPlotCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measureMergeTimeStepDataCalendarTime;
    tarch::timing::Measurement _measureMergeTimeStepDataDropFaceDataCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementAndReinitialisationCalendarTime;


  public:
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositoryArrayStack(
      peano::geometry::Geometry&                   geometry,
      int                                          maximumSizeOfCellInOutStack,
      int                                          maximumSizeOfVertexInOutStack,
      int                                          maximumSizeOfVertexTemporaryStack
    );
    
    virtual ~RepositoryArrayStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    );
         
    virtual void terminate();
        
    virtual exahype::State& getState();
    virtual const exahype::State& getState() const;

    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);
    
    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* createEmptyCheckpoint(); 

    virtual void switchToMeshRefinement();    
    virtual void switchToGridErasing();    
    virtual void switchToFusedTimeStep();    
    virtual void switchToPlotAndFusedTimeStep();    
    virtual void switchToLimiterStatusSpreading();    
    virtual void switchToReinitialisation();    
    virtual void switchToLocalRecomputationAndTimeStepSizeComputation();    
    virtual void switchToGlobalRollback();    
    virtual void switchToNeighbourDataMerging();    
    virtual void switchToSolutionUpdateAndTimeStepSizeComputation();    
    virtual void switchToTimeStepSizeComputation();    
    virtual void switchToPrediction();    
    virtual void switchToPredictionAndPlot();    
    virtual void switchToFinaliseMeshRefinementAndTimeStepSizeComputation();    
    virtual void switchToMergeTimeStepData();    
    virtual void switchToMergeTimeStepDataDropFaceData();    
    virtual void switchToFinaliseMeshRefinementAndReinitialisation();    

    virtual bool isActiveAdapterMeshRefinement() const;
    virtual bool isActiveAdapterGridErasing() const;
    virtual bool isActiveAdapterFusedTimeStep() const;
    virtual bool isActiveAdapterPlotAndFusedTimeStep() const;
    virtual bool isActiveAdapterLimiterStatusSpreading() const;
    virtual bool isActiveAdapterReinitialisation() const;
    virtual bool isActiveAdapterLocalRecomputationAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterGlobalRollback() const;
    virtual bool isActiveAdapterNeighbourDataMerging() const;
    virtual bool isActiveAdapterSolutionUpdateAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterTimeStepSizeComputation() const;
    virtual bool isActiveAdapterPrediction() const;
    virtual bool isActiveAdapterPredictionAndPlot() const;
    virtual bool isActiveAdapterFinaliseMeshRefinementAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterMergeTimeStepData() const;
    virtual bool isActiveAdapterMergeTimeStepDataDropFaceData() const;
    virtual bool isActiveAdapterFinaliseMeshRefinementAndReinitialisation() const;

     
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
