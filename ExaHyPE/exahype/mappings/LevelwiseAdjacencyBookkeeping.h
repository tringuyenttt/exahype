/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 */
#ifndef EXAHYPE_MAPPINGS_LevelwiseAdjacencyBookkeeping_H_
#define EXAHYPE_MAPPINGS_LevelwiseAdjacencyBookkeeping_H_


#include "tarch/logging/Log.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "exahype/Vertex.h"
#include "exahype/Cell.h"
#include "exahype/State.h"


namespace exahype {
      namespace mappings {
        class LevelwiseAdjacencyBookkeeping;
      } 
}


/**
 * This is an adapter which writes cell-based heap indices into the adjacency
 * maps of each vertex. This is done completely levelwisely.
 *
 * \note Hanging nodes are completely ignored and the indices initialised
 * with invalid indices everytime a hanging node is created.
 * Consistency checks must thus not consider cells which
 * have an adjacent hanging node.
 *
 * Codes which want to use this mapping need to ensure that no
 * hanging nodes appear on the boundary of the domain.
 * Refinement must be employed such that this is prevented.
 *
 * CellDescriptionsIndex   Name of the index used for the cell indices within the vertex and 
 *          the cell
 *
 * @author Tobias Weinzierl, Dominic Etienne Charrier
 * @version $Revision: 1.1 $
 */
class exahype::mappings::LevelwiseAdjacencyBookkeeping {
  public:
    peano::MappingSpecification   touchVertexLastTimeSpecification(int level) const;
    peano::MappingSpecification   touchVertexFirstTimeSpecification(int level) const;
    peano::MappingSpecification   enterCellSpecification(int level) const;
    peano::MappingSpecification   leaveCellSpecification(int level) const;
    peano::MappingSpecification   ascendSpecification(int level) const;
    peano::MappingSpecification   descendSpecification(int level) const;
    peano::CommunicationSpecification   communicationSpecification() const;

    LevelwiseAdjacencyBookkeeping();

    #if defined(SharedMemoryParallelisation)
    LevelwiseAdjacencyBookkeeping(const LevelwiseAdjacencyBookkeeping& masterThread);
    #endif

    virtual ~LevelwiseAdjacencyBookkeeping();
  
    #if defined(SharedMemoryParallelisation)
    void mergeWithWorkerThread(const LevelwiseAdjacencyBookkeeping& workerThread);
    #endif

    /**
     * Initialises the adjacency map
     * of the fine grid vertex with an invalid index.
     */
    void createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );

    /**
     * Initialises the adjacency map
     * of the fine  grid vertex with a boundary adjacency index.
     */
    void createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );

    /**
     * Initialises the adjacency map with an invalid index.
     */
    void createHangingVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );

    /**
     * Initialises the fine grid cell's heap index as invalid index.
     */
    void createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );

    /**
     * Updates the adjacency maps of all surrounding fine grid
     * vertices with the fine grid cells cell description.
     */
    void enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );

    #ifdef Parallel
    /**
     * Updates the adjacency maps of the fine grid vertices
     * at a remote boundary.
     */
    void mergeWithNeighbour(
      exahype::Vertex&  vertex,
      const exahype::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    /**
     * Updates the adjacency maps of the fine grid vertices
     * at a remote boundary. Only invoked if another mapping in
     * the same adapter requires worker->master communication.
     */
    void mergeWithMaster(
      const exahype::Cell&                     workerGridCell,
      exahype::Vertex * const                  workerGridVertices,
      const peano::grid::VertexEnumerator&     workerEnumerator,
      exahype::Cell&                           fineGridCell,
      exahype::Vertex * const                  fineGridVertices,
      const peano::grid::VertexEnumerator&     fineGridVerticesEnumerator,
      exahype::Vertex * const                  coarseGridVertices,
      const peano::grid::VertexEnumerator&     coarseGridVerticesEnumerator,
      exahype::Cell&                           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell,
      int                                      worker,
      const exahype::State&                    workerState,
      exahype::State&                          masterState
    );

    /**
     * Updates the adjacency map of a fine grid vertex
     * at a remote boundary. Only invoked if another mapping in
     * the same adapter requires master->worker communication.
     */
    void mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    /**
     * If we are on a new worker, write
     * invalid adjacency indices into the vertex'
     * adjacency map.
     */
    void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex&  localVertex,
      const exahype::Vertex&  masterOrWorkerVertex,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                       level
    );

    /**
     * If we are on a new worker, write an
     * invalid index into the cell's
     * heap index.
     */
    void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell&  localCell,
      const exahype::Cell&  masterOrWorkerCell,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                       level
    );

    void prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
    );

    bool prepareSendToWorker(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker
    );

    void prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
    );


    void receiveDataFromMaster(
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
    );


    void mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
    );


    #endif
    void destroyHangingVertex(
        const exahype::Vertex&   fineGridVertex,
        const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
        const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
        exahype::Vertex * const  coarseGridVertices,
        const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
        exahype::Cell&           coarseGridCell,
        const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );

    void destroyVertex(
        const exahype::Vertex&   fineGridVertex,
        const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
        const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
        exahype::Vertex * const  coarseGridVertices,
        const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
        exahype::Cell&           coarseGridCell,
        const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void destroyCell(
        const exahype::Cell&           fineGridCell,
        exahype::Vertex * const        fineGridVertices,
        const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
        exahype::Vertex * const        coarseGridVertices,
        const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
        exahype::Cell&                 coarseGridCell,
        const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void leaveCell(
      exahype::Cell&                          fineGridCell,
      exahype::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
    );


    void beginIteration(
      exahype::State&  solverState
    );


    void endIteration(
      exahype::State&  solverState
    );

    void descend(
      exahype::Cell * const          fineGridCells,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell
    );


    void ascend(
      exahype::Cell * const    fineGridCells,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell
    );    
};


#endif
