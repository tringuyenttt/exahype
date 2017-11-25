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
 #include "exahype/adapters/MeshRefinementAndPlotGrid.h" 
 #include "exahype/adapters/FinaliseMeshRefinement.h" 
 #include "exahype/adapters/FinaliseMeshRefinementOrLocalRollback.h" 
 #include "exahype/adapters/FusedTimeStep.h" 
 #include "exahype/adapters/PlotAndFusedTimeStep.h" 
 #include "exahype/adapters/PredictionRerun.h" 
 #include "exahype/adapters/LimiterStatusSpreading.h" 
 #include "exahype/adapters/PredictionOrLocalRecomputation.h" 
 #include "exahype/adapters/GlobalRollback.h" 
 #include "exahype/adapters/Merging.h" 
 #include "exahype/adapters/SolutionUpdate.h" 
 #include "exahype/adapters/Prediction.h" 
 #include "exahype/adapters/PredictionAndPlot.h" 



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
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MeshRefinementAndPlotGrid> _gridWithMeshRefinementAndPlotGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinement> _gridWithFinaliseMeshRefinement;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FinaliseMeshRefinementOrLocalRollback> _gridWithFinaliseMeshRefinementOrLocalRollback;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::FusedTimeStep> _gridWithFusedTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAndFusedTimeStep> _gridWithPlotAndFusedTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionRerun> _gridWithPredictionRerun;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusSpreading> _gridWithLimiterStatusSpreading;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionOrLocalRecomputation> _gridWithPredictionOrLocalRecomputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GlobalRollback> _gridWithGlobalRollback;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Merging> _gridWithMerging;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionUpdate> _gridWithSolutionUpdate;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Prediction> _gridWithPrediction;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndPlot> _gridWithPredictionAndPlot;

  
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureMeshRefinementCPUTime;
    tarch::timing::Measurement _measureMeshRefinementAndPlotGridCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementCPUTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementOrLocalRollbackCPUTime;
    tarch::timing::Measurement _measureFusedTimeStepCPUTime;
    tarch::timing::Measurement _measurePlotAndFusedTimeStepCPUTime;
    tarch::timing::Measurement _measurePredictionRerunCPUTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCPUTime;
    tarch::timing::Measurement _measurePredictionOrLocalRecomputationCPUTime;
    tarch::timing::Measurement _measureGlobalRollbackCPUTime;
    tarch::timing::Measurement _measureMergingCPUTime;
    tarch::timing::Measurement _measureSolutionUpdateCPUTime;
    tarch::timing::Measurement _measurePredictionCPUTime;
    tarch::timing::Measurement _measurePredictionAndPlotCPUTime;

    tarch::timing::Measurement _measureMeshRefinementCalendarTime;
    tarch::timing::Measurement _measureMeshRefinementAndPlotGridCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementCalendarTime;
    tarch::timing::Measurement _measureFinaliseMeshRefinementOrLocalRollbackCalendarTime;
    tarch::timing::Measurement _measureFusedTimeStepCalendarTime;
    tarch::timing::Measurement _measurePlotAndFusedTimeStepCalendarTime;
    tarch::timing::Measurement _measurePredictionRerunCalendarTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCalendarTime;
    tarch::timing::Measurement _measurePredictionOrLocalRecomputationCalendarTime;
    tarch::timing::Measurement _measureGlobalRollbackCalendarTime;
    tarch::timing::Measurement _measureMergingCalendarTime;
    tarch::timing::Measurement _measureSolutionUpdateCalendarTime;
    tarch::timing::Measurement _measurePredictionCalendarTime;
    tarch::timing::Measurement _measurePredictionAndPlotCalendarTime;


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
    virtual void switchToMeshRefinementAndPlotGrid();    
    virtual void switchToFinaliseMeshRefinement();    
    virtual void switchToFinaliseMeshRefinementOrLocalRollback();    
    virtual void switchToFusedTimeStep();    
    virtual void switchToPlotAndFusedTimeStep();    
    virtual void switchToPredictionRerun();    
    virtual void switchToLimiterStatusSpreading();    
    virtual void switchToPredictionOrLocalRecomputation();    
    virtual void switchToGlobalRollback();    
    virtual void switchToMerging();    
    virtual void switchToSolutionUpdate();    
    virtual void switchToPrediction();    
    virtual void switchToPredictionAndPlot();    

    virtual bool isActiveAdapterMeshRefinement() const;
    virtual bool isActiveAdapterMeshRefinementAndPlotGrid() const;
    virtual bool isActiveAdapterFinaliseMeshRefinement() const;
    virtual bool isActiveAdapterFinaliseMeshRefinementOrLocalRollback() const;
    virtual bool isActiveAdapterFusedTimeStep() const;
    virtual bool isActiveAdapterPlotAndFusedTimeStep() const;
    virtual bool isActiveAdapterPredictionRerun() const;
    virtual bool isActiveAdapterLimiterStatusSpreading() const;
    virtual bool isActiveAdapterPredictionOrLocalRecomputation() const;
    virtual bool isActiveAdapterGlobalRollback() const;
    virtual bool isActiveAdapterMerging() const;
    virtual bool isActiveAdapterSolutionUpdate() const;
    virtual bool isActiveAdapterPrediction() const;
    virtual bool isActiveAdapterPredictionAndPlot() const;

     
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
