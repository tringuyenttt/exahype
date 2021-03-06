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
 **/
 
#ifndef _EXAHYPE_CELL_H_
#define _EXAHYPE_CELL_H_

#include "exahype/State.h"

#include <deque>

#include "peano/grid/Cell.h"
#include "peano/grid/VertexEnumerator.h"

#include "peano/utils/Globals.h"

#include "exahype/records/Cell.h"

#include "exahype/records/ADERDGCellDescription.h"
#include "exahype/records/FiniteVolumesCellDescription.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

namespace exahype {
  class Cell;

  namespace solvers {
    class ADERDGSolver;
    class FiniteVolumesSolver;
  }
}

/**
 * @todo Dominic We should add some proper descriptions here one day.
 */
class exahype::Cell : public peano::grid::Cell<exahype::records::Cell> {
 private:
  typedef class peano::grid::Cell<exahype::records::Cell> Base;

  static tarch::logging::Log _log;

  #ifdef Parallel
  /**
   * Heap index for tempoarily storing metadata between
   * the calls of Mapping::receiveDataFromMaster and Mapping::mergeWithWorker.
   */
  static int ReceivedMetadataHeapIndex;

  /**
   * Heap index for temporarily storing metadata between
   * Mapping::receiveDataFromMaster and Mapping::mergeWithWorker.
   */
  static std::deque<int> ReceivedDataHeapIndices;
  #endif

 public:
  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  Cell();

  /**
   * This constructor should not set any attributes. It is used by the
   * traversal algorithm whenever it allocates an array whose elements
   * will be overwritten later anyway.
   */
  Cell(const Base::DoNotCallStandardConstructor&);

  /**
   * Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it. It is kind of a copy constructor that converts an object which
   * comprises solely persistent attributes into a full attribute. This very
   * functionality is implemented within the super type, i.e. this constructor
   * has to invoke the correponsing super type's constructor and not the super
   * type standard constructor.
   */
  Cell(const Base::PersistentCell& argument);

  /**
   * Validates that all incoming neighbour merges have been
   * performed.
   *
   * <h2> ADER-DG and Limiting-ADERDG </h2>
   *
   * Descendants might border to uninitialised Peano cells.
   * A merge is thus not performed here.
   */
  static void validateThatAllNeighbourMergesHaveBeenPerformed(
      const int cellDescriptionsIndex,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  /**
   * Here we reset helper variables that play a role in
   * the neighbour merge methods.
   * These are the cell description attributes
   * riemannSolvePerformed[DIMENSIONS_TIMES_TWO], and
   * faceDataExchangeCounter[DIMENSIONS_TIMES_TWO].
   *
   * <h2>Shared Memory</h2>
   * The flag riemannSolvePerformed
   * indicates for every thread that touches a
   * face of a cell description if a Riemann Solve
   * was already performed for this face.
   *
   * <h2>MPI</h2>
   * This method resets Face data exchange counters:
   * To this end, we count the listings of a remote rank on each
   * of the faces surrounding a cell description.
   *
   * We perform the following actions depending on the counter value:
   * 0 - no connection: no send. Set to unreachable value.
   * 2^{d-2} - full face connection where cell is inside but half of face vertices are outside:
   * send at time of 2^{d-2}-th touch of face.
   * 4^{d-2} - full face connection where cell is inside and face vertices are all inside:
   * send at time of 2^{d-2}-th touch of face.
   */
  static void resetNeighbourMergeFlags(
      const int cellDescriptionsIndex);

  /**
   * Here we reset helper variables that play a role in
   * the neighbour merge methods.
   * These are the cell description attributes
   * riemannSolvePerformed[DIMENSIONS_TIMES_TWO], and
   * faceDataExchangeCounter[DIMENSIONS_TIMES_TWO].
   *
   * <h2>Shared Memory</h2>
   * The flag riemannSolvePerformed
   * indicates for every thread that touches a
   * face of a cell description if a Riemann Solve
   * was already performed for this face.
   *
   * <h2>MPI</h2>
   * This method resets Face data exchange counters:
   * To this end, we count the listings of a remote rank on each
   * of the faces surrounding a cell description.
   *
   * We perform the following actions depending on the counter value:
   * 0 - no connection: no send. Set to unreachable value.
   * 2^{d-2} - full face connection where cell is inside but half of face vertices are outside:
   * send at time of 2^{d-2}-th touch of face.
   * 4^{d-2} - full face connection where cell is inside and face vertices are all inside:
   * send at time of 2^{d-2}-th touch of face.
   */
  static void resetFaceDataExchangeCounters(
      const int cellDescriptionsIndex,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  /**
   * Determine inside and outside faces of a cell.
   * A face is considered inside if at least
   * one of its vertices is inside.
   *
   * <h2>Issues with AMR</h2>
   * Use this function only for cells
   * on the coarsest mesh level per solver.
   *
   * If the computational domain
   * is anisotropic and we perform adaptive refinement.
   * Then, fine grid cells might suddenly get neighbours and thus
   * their vertices are not outside anymore.
   *
   * For newly introduced fine mesh cells, then prolongate the coarse grid
   * information down.
   */
  static std::bitset<DIMENSIONS_TIMES_TWO> determineInsideAndOutsideFaces(
        const exahype::Vertex* const verticesAroundCell,
        const peano::grid::VertexEnumerator& verticesEnumerator);

  /**
   * Computes the barycentre of a face of a cell.
   *
   * \param[in] cellDescription cell description where we get the geometry information from
   * \param[in] direction       coordinate direction the face's outward directed normal vector is aligned with (0:x, 1:y, 2: z)
   * \param[in] orientation     orientation of the face's outward directed normal vector (0: -1, 1: +1)
   */
  static tarch::la::Vector<DIMENSIONS,double> computeFaceBarycentre(
      const tarch::la::Vector<DIMENSIONS,double>& cellOffset,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int direction, const int orientation);

  /**
   * Returns true if the cell corresponding
   * to the vertices \p verticesAroundCell
   * is neighbour to a remote rank
   * via one of the faces.
   * Only inside faces are checked, i.e. faces where
   * at least one vertex is inside.
   * Boundary and outside vertices are ignored.
   */
  static bool isAdjacentToRemoteRankAtInsideFace(
      exahype::Vertex* const               verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);

  /**
   * Returns meta data describing the surrounding cell descriptions. The
   * routine is notably used by the automated adapters to derive adjacency
   * information on the cell level.
   */
  int getCellDescriptionsIndex() const;

  /**
   * TODO(Dominic): Add docu.
   */
  void setCellDescriptionsIndex(int cellDescriptionsIndex);

  /**
   * Returns the number of ADERDGCellDescriptions associated
   * with this cell.
   */
  int getNumberOfADERDGCellDescriptions() const;

  /**
   * Returns the number of FiniteVolumesCellDescriptions associated
   * with this cell.
   */
  int getNumberOfFiniteVolumeCellDescriptions() const;

  /**
   * Each cell points to a series of cell descriptions. The array holding the
   * series has to be stored on the heap and, consequently, initialised
   * properly. This is done by create() while destroy() cleans up. Please note
   * that you have to invoke create() once before you do anything with the cell
   * at all. You should destroy() in return in the very end.
   *
   * The operation shows that each cell in the tree can theoretically hold a
   * solver though only few do.
   *
   * This operation is used by addNewCellDescription().
   */
  void setupMetaData();

  /**
   * @see setupMetaData()
   */
  void shutdownMetaData();

  /**
   * \return true if no cell descriptions
   * are allocated for this cell.
   */
  bool isEmpty() const;

  /**
   * todo docu
   *
   * setupMetaData() is called if cell hasn't been properly initialised before.
   *
   * \note Operation is thread-safe.
   */
  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::ADERDGCellDescription::Type cellType,
      const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
      const int level,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const tarch::la::Vector<DIMENSIONS, double>& cellOffset);

  /**
    * TODO(Dominic): Docu.
    *
    *  \note Operation is thread-safe.
    */
  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::FiniteVolumesCellDescription::Type cellType,
      const exahype::records::FiniteVolumesCellDescription::RefinementEvent refinementEvent,
      const int level,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const tarch::la::Vector<DIMENSIONS, double>& cellOffset);

  /**
   * @return if this cell is initialised.
   *
   * @developers:
   * Note that it is not simply sufficient to
   * check if the heap index equals
   * multiscalelinkedcell::HangingVertexBookkeeper.
   *
   * @TODO bug
   */
  bool isInitialised() const;

  #ifdef Parallel
  /**
   * \return true if the cell is inside and the
   * cell size does belong to a grid level that is
   * occupied by a solver.
   *
   * \note isInside() is not necessarily consistent on receivedCell from receiveDataFromMaster and
   * on localCell in mergeWithWorker and thus ignored here.
   */
  bool hasToCommunicate(
      const tarch::la::Vector<DIMENSIONS,double>& cellSize ) const;

  /**
   * Count the listings of remote ranks sharing an inside vertex
   * adjacent to the face \p faceIndex of a cell with this rank.
   * In case all vertices adjacent to the face are inside of the domain,
   * this value is either 0 or 2^{d-1}.
   *
   * If we count 2^{d-1} listings, we directly know that this rank
   * shares a whole face with a remote rank.
   * If we count 0 listings, we do not have
   * a remote rank adjacent to this face.
   *
   * We know from the result of this function how
   * many vertices will try to exchange neighbour information
   * at this face.
   *
   * <h2>Boundary vertices</h2>
   * Boundary vertices are ignored by the counting since they are
   * skipped in all merging and sending routines;
   * see method exahype::Vertex::hasToCommunicate.
   * This introduces further values for
   * the counted listings:
   *
   * 2^{d-2} - two vertices belong to the boundary.
   * 2^{d-3} - (only in 3d) three vertices belong to the boundary. This might
   *           happen in a corner of the domain if the bounding
   *           box is not scaled.
   *
   * @developers:
   * TODO(Dominic): We currently check for uniqueness of the
   * remote rank. This might however not be necessary.
   */
  static int countListingsOfRemoteRankAtInsideFace(
      const int                            faceIndex,
      exahype::Vertex* const               verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);

  // MASTER->WORKER

  /*!
   * Send data such global solver and plotter
   * time step data down to a worker.
   */
  static void broadcastGlobalDataToWorker(
      const int                                   worker,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const int                                   level);

  /*!
   * Merge with global data such as global solver and plotter
   * time step data sent down from the master.
   */
  static void mergeWithGlobalDataFromMaster(
      const int                                   master,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const int                                   level);

  /*! Broadcast metadata down to the worker.
   * This is done per cell.
   */
  void broadcastMetadataToWorkerPerCell(
      const int                                   worker,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level);

  /*! Receive the metadata broadcasted down from the master.
   * This is done per cell.
   *
   * \note must be called on the "receivedCell" in
   * "Mapping::receiveDataFromMaster".
   */
  void receiveMetadataFromMasterPerCell(
      const int                                   master,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level);

  /*! Receives metadata from the master and merges it with all
   * solvers registered on the cell.
   *
   * \note must be called on the "localCell" in
   * "Mapping::mergeWithWorker".
   */
  void mergeWithMetadataFromMasterPerCell(
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const exahype::State::AlgorithmSection&     section);

  /*! Broadcast data per cell to a worker.
   *
   * Loop over all solvers and send data down
   * to the worker if the respective solver
   * has registered a patch at this cell which
   * requires such a send.
   * In all other cases, send a zero-length
   * message.
   *
   * \param[in] worker the worker rank.
   * \param[in] cellCentre centre of this cell.
   * \param[in] cellSize   size of this cell.
   * \param[in] level      grid level this cell is residing at.
   */
  void broadcastDataToWorkerPerCell(
    const int                                   worker,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize,
    const int                                   level) const;

  /*! Receive cell-wise heap data from the master.
   *
   * \note Must be called on the "receivedCell" in "Mapping::receiveDataFromMaster".
   *
   * \param[in] cellCentrie centre of the received cell.
   * \param[in] cellSize   size of this cell.
   * \param[in] level grid level the received cell resides at.
   *
   */
  void receiveDataFromMasterPerCell(
      const int                                   master,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level);

  /*! Merge received heap data with the local cell.
   *
   * In Peano's two-step receive-from-master process,
   * we first receive data in "Mapping::receiveDataFromMaster".
   * We store received heap data in the fields ReceivedDataHeapIndex and
   * ReceivedHeapDataIndices.
   *
   * In the second and last step, in "Mapping::mergeWithWorker",
   * we pick up the previously received heap data, merge
   * it with the "localCell", and then delete the heap data associated
   * with the "receivedCell".
   *
   * \note Must be called on the "localCell" in "Mapping::mergeWithWorker".
   *
   * \param[in] cellSize size of either cell.
   */
  void mergeWithMasterDataPerCell(
      const tarch::la::Vector<DIMENSIONS,double>& cellSize );

  // WORKER->MASTER

  /**
   * Reduce Metadata from a worker
   * to the master.
   */
  void reduceMetadataToMasterPerCell(
      const int                                   master,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level) const;

  /**
   * Receives metadata from a worker and merges it with all
   * solvers registered on the cell.
   *
   * \return if vertical (master-worker) exchange
   * of face data is required during the time stepping iterations
   * for any of the registered solvers.
   */
  bool mergeWithMetadataFromWorkerPerCell(
      const int                                   workerRank,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level,
      const exahype::State::AlgorithmSection&     section) const;


  /**
   * Reduce metadata and face data to the master
   * In contrast to the two-step broadcast merging
   * process, reductions are a single-step process.
   *
   * \param[in] master     the master rank.
   * \param[in] cellCentre centre of this cell.
   * \param[in] cellSize   size of this cell.
   * \param[in] level      grid level this cell is residing at.
   */
  void reduceDataToMasterPerCell(
      const int                                   master,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level) const;

  /**
   * Merge metadata and face data from the worker.
   * In contrast to the two-step broadcast merging
   * process, reductions are a single-step process.
   *
   * \param[in] worker     the worker rank.
   * \param[in] cellCentre centre of this cell.
   * \param[in] cellSize   size of this cell.
   * \param[in] level      grid level this cell is residing at.
   */
  void mergeWithDataFromWorkerPerCell(
      const int                                   worker,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>& cellSize,
      const int                                   level) const;

  // global

  /*!
   * Send data such global solver
   * time step data up to the master.
   */
  static void reduceGlobalDataToMaster(
    const int                                   master,
    const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
    const int                                   level);

  /*!
   * Merge with global data such as global solver
   * time step data sent up from the master.
   */
  static void mergeWithGlobalDataFromWorker(
      const int                                   worker,
      const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
      const int                                   level);
  #endif
};

#endif
