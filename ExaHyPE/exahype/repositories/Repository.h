// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPEREPOSITORIES__REPOSITORY_H_ 
#define _EXAHYPEREPOSITORIES__REPOSITORY_H_


#include <string>
#include <vector>

#include "tarch/logging/Log.h"
#include "tarch/timing/Measurement.h"

#include "exahype/State.h"

#include "peano/grid/Checkpoint.h"


namespace exahype {
      namespace repositories {
        class Repository;
      } 
}



/**
 * Interface of all repositories
 */
class exahype::repositories::Repository {
  public:
    virtual ~Repository() {};
    
    /**
     * Iterate with current active event handle.
     *
     * <h2> Stop exchange of boundary vertices </h2>
     *
     * If you can make your algorithm run without boundary vertex exchange for
     * a couple of times, you typically observe a significant speedup. Exchange
     * via a heap in return can often be hidden in the background easier. If you
     * unset the argument exchangeBoundaryVertices, this has two implications:
     *
     * - All vertices are passed through perpareSendToNeighbour and merge to
     *   allow you to do heap data exchange or boundary modifications, e.g.
     *   However, the vertices are not physically handled over the MPI but
     *   discarded instead.
     * - If you switch to false/true this always has implications on the next
     *   grid traversal, as boundary data is exchanged similar to Jacobi
     *   sweeps. If you send away data and afterward switch the flag to false,
     *   there is still data available in the next sweep that has been received
     *   from the iteration before. You do not send away data anymore, but you
     *   receive.
     *
     * <h2>iterate on MPI workers</h2>
     *
     * If iterate is called on MPI workers, i.e. not on the global rank, then 
     * the two parameters are actually ignored. Instead, the iterate relies on 
     * the existing _repositoryState that has been received by 
     * continueToIterate() before. 
     *
     * @param numberOfIterations        How often shall the repository run through a particular adapter
     * @param exchangeBoundaryVertices  Does Peano have to exchange boundary vertices at all?
     */
    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true) = 0;

    virtual exahype::State& getState() = 0;
    virtual const exahype::State& getState() const = 0;

    /**
     * Switch to another event handle.
     */
    virtual void switchToMeshRefinement() = 0;    
    virtual void switchToMeshRefinementAndPlotGrid() = 0;    
    virtual void switchToFinaliseMeshRefinement() = 0;    
    virtual void switchToFinaliseMeshRefinementOrLocalRollback() = 0;    
    virtual void switchToFusedTimeStep() = 0;    
    virtual void switchToPredictionRerun() = 0;    
    virtual void switchToBroadcastGlobalDataAndDropNeighbourMessages() = 0;    
    virtual void switchToLimiterStatusSpreading() = 0;    
    virtual void switchToPredictionOrLocalRecomputation() = 0;    
    virtual void switchToGlobalRollback() = 0;    
    virtual void switchToBroadcastGlobalDataAndMergeNeighbourMessages() = 0;    
    virtual void switchToSolutionUpdate() = 0;    
    virtual void switchToPrediction() = 0;    

    virtual bool isActiveAdapterMeshRefinement() const = 0;
    virtual bool isActiveAdapterMeshRefinementAndPlotGrid() const = 0;
    virtual bool isActiveAdapterFinaliseMeshRefinement() const = 0;
    virtual bool isActiveAdapterFinaliseMeshRefinementOrLocalRollback() const = 0;
    virtual bool isActiveAdapterFusedTimeStep() const = 0;
    virtual bool isActiveAdapterPredictionRerun() const = 0;
    virtual bool isActiveAdapterBroadcastGlobalDataAndDropNeighbourMessages() const = 0;
    virtual bool isActiveAdapterLimiterStatusSpreading() const = 0;
    virtual bool isActiveAdapterPredictionOrLocalRecomputation() const = 0;
    virtual bool isActiveAdapterGlobalRollback() const = 0;
    virtual bool isActiveAdapterBroadcastGlobalDataAndMergeNeighbourMessages() const = 0;
    virtual bool isActiveAdapterSolutionUpdate() const = 0;
    virtual bool isActiveAdapterPrediction() const = 0;


    /**
     * Give Some Statistics
     *
     * This operation gives you a table which tells you for each adapter how 
     * much time was spent in it. The result is written to the info log device. 
     *
     * @param logAllAdapters If it is set to true then all adapters
     *          are shown. If you hand in false, then only those adapters that
     *          have been in use are displayed.
     */
    virtual void logIterationStatistics(bool logAllAdapters) const = 0;
        
    virtual void clearIterationStatistics() = 0;

    /**
     * Create a checkpoint.
     *
     * See createEmptyCheckpoint() before.
     */
    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex,exahype::Cell> * const checkpoint) = 0; 
    
    /**
     * Load a checkpoint
     * 
     * Does neither modify the checkpoint nor does it delete it. If you want to 
     * load a file from a checkpoint, see createEmptyCheckpoint() before.
     */
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex,exahype::Cell> const * const checkpoint ) = 0;
    
    /**
     * Create empty Checkpoint
     *
     * If you wanna read a checkpoint, implement the following four steps:
     * - Call createEmptyCheckpoint() on the repository. You receive a pointer 
     *   to a new checkpoint object. If you don't use this operation, your code 
     *   won't work in parallel and is not grid-independent.
     * - Invoke readFromFile() on the checkpoint object.
     * - Call readCheckpoint() on the repository and pass it your new checkpoint 
     *   object.
     * - Destroy the checkpoint object on the heap. 
     *
     * If you wanna write a checkpoint, implement the following four steps:
     * - Call createEmptyCheckpoint() on the repository. You receive a pointer 
     *   to a new checkpoint object. If you don't use this operation, your code 
     *   won't work in parallel and is not grid-independent.
     * - Call writeCheckpoint() on the repository and pass it your new checkpoint 
     *   object.
     * - Invoke writeToFile() on the checkpoint object.
     * - Destroy the checkpoint object on the heap. 
     */
    virtual peano::grid::Checkpoint<exahype::Vertex,exahype::Cell>* createEmptyCheckpoint() = 0;
    
    /**
     * Restart the repository with a different setup.
     *
     * This operation is only used by the parallel code. It itself derives from 
     * the new master node the new state and the adjacency information. Also 
     * the vertices with the adjacency information are taken from the master 
     * node, so the only interesting thing is how to traverse the adjacent 
     * elements and how the curve runs through the root element.
     */
    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    ) = 0;     
    
    /**
     * Counterpart of restart(). However, terminate() also is to be called on  
     * the global master.
     */
    virtual void terminate() = 0;
    
    #ifdef Parallel
    enum ContinueCommand {
      Continue,
      Terminate,
      RunGlobalStep
    };
    
    /**
     * Shall a worker in the parallel cluster continue to iterate?
     *
     * This operation may be invoked on a worker node only, i.e. you are not 
     * allowed to trigger it on global rank 0. It waits for a wake-up call from 
     * the master node and then tells you whether to continue your work or not. 
     * If the result is false, you might do some additional iterations (plotting 
     * stuff or writing some statistics, e.g.), but then you should call
     * terminate() on the node and ask the node pool for a new job. If you 
     * invoke additional iterates() after this operation has returned false and 
     * and the terminate(), these iterations won't trigger any communication 
     * anymore.
     */
    virtual ContinueCommand continueToIterate() = 0;

    /**
     * Run one global step on all mpi ranks besides rank 0
     * 
     * This operation sends a marker to all nodes, i.e. 
     * both idle and working nodes, and calls their runGlobalStep() routine 
     * within the parallel runner. Afterwards, all idle nodes again register as 
     * idle on the node pool, all other nodes continue to run Peano. Should be 
     * used with care, as it might be expensive on massively parallel systems.  
     *
     * Please note that all worker runners invoke their runGlobalStep() as a 
     * reaction to this function call. The global master, i.e. rank 0, however
     * does not run this operation. If you want to make the master trigger the 
     * function as well, you have to invoke it manually. Typically, the global 
     * master code thus comprises two steps: call runGlobalStep() on the 
     * repository (which in turn invokes runGlobalStep() on all workers) and 
     * immediately afterwards call runGlobalStep() locally within 
     * runAsMaster().
     *
     * !!! Dynamic load balancing
     *
     * This operation should not be used if you have dynamic load balancing 
     * enabled. If any rank combination currently is involved in a join or 
     * fork, the worker ranks do not receive the global step instruction. If 
     * you want to use this operation, ensure no load balancing is happening:
     * either through your application-specific knowledge or by switching it 
     * manually off. Please note that switching load balancing off might mean 
     * that you still need up to three time steps to finish all running load 
     * balancing activities.
     */
    virtual void runGlobalStep() = 0;
    #endif
    
    
    /**
     * Set maximum memory footprint spent on temporary data
     *
     * This value by default is infinity and, thus, no manual restriction on 
     * the maximum memory footprint spent on temporary data is posed. And, 
     * theoretically, the more memory you allow Peano to spend on temporary 
     * data the better the parallel shared memory scalability, as the 
     * concurrency level raises if more data is held temporary. However, due 
     * to NUMA effects, sometimes codes perform better if you restrict this 
     * value. Also, you might wanna run into swapping if this value is not 
     * set manually.
     */    
     virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value) = 0;    
};


#endif
