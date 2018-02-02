#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexSTDStack.h"

 #include "exahype/adapters/MeshRefinement.h" 
 #include "exahype/adapters/MeshRefinementAndPlotGrid.h" 
 #include "exahype/adapters/FinaliseMeshRefinement.h" 
 #include "exahype/adapters/FinaliseMeshRefinementOrLocalRollback.h" 
 #include "exahype/adapters/FusedTimeStep.h" 
 #include "exahype/adapters/PredictionRerun.h" 
 #include "exahype/adapters/BroadcastGlobalDataAndDropNeighbourMessages.h" 
 #include "exahype/adapters/LimiterStatusSpreading.h" 
 #include "exahype/adapters/PredictionOrLocalRecomputation.h" 
 #include "exahype/adapters/GlobalRollback.h" 
 #include "exahype/adapters/BroadcastGlobalDataAndMergeNeighbourMessages.h" 
 #include "exahype/adapters/SolutionUpdate.h" 
 #include "exahype/adapters/Prediction.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::PredictionOrLocalRecomputation>;
  }
}

#include "peano/grid/Grid.cpph"
