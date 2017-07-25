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
 
#ifndef _EXAHYPE_VERTEX_H_
#define _EXAHYPE_VERTEX_H_

#include "exahype/records/Vertex.h"
#include "peano/grid/Vertex.h"
#include "peano/grid/VertexEnumerator.h"
#include "peano/utils/Globals.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

namespace exahype {
/**
 * We abuse this heap to send and receive metadata from one MPI rank to the other.
 * We never actually store data on this heap.
 * TODO(Dominic): Change to RLEIntegerHeap that compresses data.
 */
typedef peano::heap::PlainIntegerHeap  MetadataHeap; // TODO(Dominic): Migrate to Vertex.

class Vertex;

/**
 * Forward declaration
 */
class VertexOperations;
}

/**
 * A grid vertex.
 *
 * Peano realises the neighbour communication between
 * different MPI ranks via exchanging and merging vertices.
 * It further offers routines in the mappings to prepare
 * vertices before the data exchange and routines to merge
 * the exchanged vertices. In these two routines, we plugin the
 * sending and receiving of heap data.
 * A fair share of the required heap data exchange functionality can
 * thus be found in this class.
 */
class exahype::Vertex : public peano::grid::Vertex<exahype::records::Vertex> {
 private:
  typedef class peano::grid::Vertex<exahype::records::Vertex> Base;

  friend class VertexOperations;

  static tarch::logging::Log _log;
 public:

  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  Vertex();

  /**
   * This constructor should not set any attributes. It is used by the
   * traversal algorithm whenever it allocates an array whose elements
   * will be overwritten later anyway.
   */
  Vertex(const Base::DoNotCallStandardConstructor&);

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
  Vertex(const Base::PersistentVertex& argument);

  /**
   * Return the cell descriptions indices of the adjacent cells.
   */
  tarch::la::Vector<TWO_POWER_D, int>& getCellDescriptionsIndex();

//  struct Face {
//    int normalDirection;
//    int faceIndex1;
//    int faceIndex2
//  } typedef Face;
//
//  static Face computeFaceIndices(
//      const tarch::la::Vector<DIMENSIONS,int>& pos1,
//      const tarch::la::Vector<DIMENSIONS,int>& pos2) {
//    assertion(tarch::la::countEqualEntries(pos1,pos2)==DIMENSIONS-1);
//    assertion(tarch::la::equalsReturnIndex(pos1, pos2)<DIMENSIONS);
//    Face face;
//
//    face.normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
//    assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
//    face.faceIndex1 = 2 * face.normalDirection +
//        (pos2(face.normalDirection) > pos1(face.normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
//    face.faceIndex2      = 2 * face.normalDirection +
//        (pos1(face.normalDirection) > pos2(face.normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
//
//    assertion(face.normalDirection<DIMENSIONS);
//    assertion(std::abs(face.faceIndex1-face.faceIndex2)==1);
//    return face;
//  }

  /**
   * Checks if the cell descriptions at the indices corresponding
   * to \p pos1 and \p pos2 need to be merged with each other.
   *
   * TODO(Dominic): The idea is to store purely geometry based information
   * (offset,size,riemannSolvePerfomed,..) on a separate heap.
   * That's why I have not merged the loops in this method
   * into the solvers. I need to discuss this with Tobias.
   */
  bool hasToMergeNeighbours(
        const tarch::la::Vector<DIMENSIONS,int>& pos1,
        const tarch::la::Vector<DIMENSIONS,int>& pos2) const;

  /**
   * Checks if the cell description at the indices corresponding
   * to \p pos1 and \p pos2 need to be merged with each other.
   *
   * TODO(Dominic): The idea is to store purely geometry based information
   * (offset,size,riemannSolvePerfomed,..) on a separate heap.
   * That's why I have not merged the loops in this method
   * into the solvers. I need to discuss this with Tobias.
   */
  bool hasToMergeWithBoundaryData(
        const tarch::la::Vector<DIMENSIONS,int>& pos1,
        const tarch::la::Vector<DIMENSIONS,int>& pos2) const;

  /**
   * Sets a flag on the cell descriptions at the indices corresponding
   * to \p pos1 and \p pos2 that the merge with the neighbours has been
   * performed.
   *
   * TODO(Dominic): The idea is to store purely geometry based information
   * (offset,size,riemannSolvePerfomed,..) on a separate heap.
   * That's why I have not merged the loops in this method
   * into the solvers. I need to discuss this with Tobias.
   */
  void setMergePerformed(
          const tarch::la::Vector<DIMENSIONS,int>& pos1,
          const tarch::la::Vector<DIMENSIONS,int>& pos2,
          bool state) const;


#ifdef Parallel
  /**
   * Defines an invalid metadata entry.
   */
  static const int InvalidMetadataEntry;

  /**
   * Defines the length of the metadata
   * we send out per sovler.
   */
  static const int MetadataPerSolver;

  /**
   * Encodes the metadata as integer sequence.
   *
   * The first element refers to the number of
   * ADERDGCellDescriptions associated with this cell (nADERG).
   * The next 2*nADERG elements store a pair of
   * solver number, and cell description type (encoded as int)
   * for each ADERDGCellDescription associated with this cell (description).
   *
   * The element 1+2*nADERDG refers to the number of
   * FiniteVolumesCellDescriptions associated with this cell (nFV).
   * The remaining 2*nFV elements store a pair of
   * solver number, and cell description type (encoded as int)
   * for each FiniteVolumesCellDescription associated with this cell
   * (description).
   *
   * @developers:
   * TODO(Dominic): Does it make sense to also encode the
   *                refinement event?
   * TODO(Dominic): Not directly associated with a cell. Consider
   * to move this function somewhere else.
   */
  static exahype::MetadataHeap::HeapEntries encodeMetadata(const int cellDescriptionsIndex);

  /**
   * Creates a sequence of \p InvalidMetadataEntry with length
   * exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver.
   */
  static exahype::MetadataHeap::HeapEntries createEncodedMetadataSequenceWithInvalidEntries() {
      exahype::MetadataHeap::HeapEntries encodedMetaData(
          exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver,
          exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver);
      std::fill_n(encodedMetaData.begin(),encodedMetaData.size(),InvalidMetadataEntry); // Implicit conversion.
      return encodedMetaData;
  }

  /**
   * Checks if all the entries of \p sequence are set to
   * \p InvalidMetadataEntry.
   */
  static bool isEncodedMetadataSequenceWithInvalidEntries(exahype::MetadataHeap::HeapEntries& sequence) {
     assertion(sequence.size() == exahype::solvers::RegisteredSolvers.size()*MetadataPerSolver);

     for (auto& m : sequence)
       if (m.getU()==InvalidMetadataEntry)
         return false;

     return true;
  }

  /**
   * Send metadata to rank \p toRank.
   */
  static void sendEncodedMetadata(
      const int                                   toRank,
      const int                                   cellDescriptionsIndex,
      const peano::heap::MessageType&             messageType,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level) {
    MetadataHeap::HeapEntries encodedMetadata =
        encodeMetadata(cellDescriptionsIndex);
    MetadataHeap::getInstance().sendData(
        encodedMetadata,toRank,x,level,messageType);
  }

  /**
   * Send a metadata sequence filled with InvalidMetadataEntry
   * to rank \p toRank.
   */
  static void sendEncodedMetadataSequenceWithInvalidEntries(
      const int                                   toRank,
      const peano::heap::MessageType&             messageType,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level) {
    MetadataHeap::HeapEntries encodedMetadata =
        createEncodedMetadataSequenceWithInvalidEntries();
    MetadataHeap::getInstance().sendData(
        encodedMetadata,toRank,x,level,messageType);
  }

  /**
   * Drop metadata sent by rank \p fromRank.
   */
  static void dropMetadata(
      const int                                   fromRank,
      const peano::heap::MessageType&             messageType,
      const tarch::la::Vector<DIMENSIONS,double>& x,
      const int                                   level) {
    MetadataHeap::getInstance().receiveData(
        fromRank,x,level,messageType);
  }

  /**
   * Returns if this vertex needs to send a metadata message to a remote rank \p toRank.
   *
   * We need to send a message to remote rank \p roRank if both ranks
   * share a face and the adjacent rank at position dest is the remote rank.
   *
   * It is further necessary that either holds:
   *
   * 1. The adjacent rank at position src in the vertex' adjacency information
   *    equals the rank of the MPI process that calls
   *    this function.
   *
   * 2. For the rank at position src in the vertex' adjacency information forking was triggered.
   *    Then, the domain at position src will not be owned anymore by the rank that calls
   *    this function in the next iteration. However, remote ranks still expect receiving data from it.
   *
   * @param state  The state tells us if a neighbouring rank is forking or if forking was triggered for this rank.
   * @param src    A counter of a d-dimensional for-loop referring to the the message source
   *               in the vector returned by getAdjacentRemoteRanks().
   * @param dest   A counter of a d-dimensional for-loop referring to the message destination
   *               in the vector returned by getAdjacentRemoteRanks().
   * @param toRank The rank we want to send the message to.
   *
   * @developers:
   * TODO(Dominic): Consider joins.
   * TODO(Dominic): Potentially, there is still a bug if two neighbouring ranks are
   *                forking at the same time.
   */
  bool hasToSendMetadata(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int toRank);

  /**
   * Similar to hasToSendMetadata. However ignores
   * that the dest rank in the adjacency information
   * might be a forking/joining one.
   */
  bool hasToSendMetadataIgnoreForksJoins(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int toRank);

  /**
   * Returns if this vertex needs to receive a metadata message from a remote rank \p fromRank.
   *
   * We need to receive a message from remote rank \p fromRank if both ranks
   * share a face and the adjacent rank at position src is the remote rank.
   *
   * It is further necessary that either holds:
   *
   * 1. The adjacent rank at position dest in the vertex' adjacency information
   *    equals the rank of the MPI process that calls
   *    this function.
   *
   * 2. The rank at position dest in the vertex' adjacency information is now a forking rank.
   *    Then, the domain at position dest was owned by the rank of the MPI process that calls
   *    this function in the previous iteration and remote ranks have thus sent data to it.
   *
   *
   * @param state    The state tells us if a neighbouring rank is forking or if forking was triggered for this rank.
   * @param src      A counter of a d-dimensional for-loop referring to the the message source
   *                 in the vector returned by getAdjacentRemoteRanks().
   * @param dest     A counter of a d-dimensional for-loop referring to the message destination
   *                 in the vector returned by getAdjacentRemoteRanks().
   * @param fromRank The rank we want to send the message to.
   *
   * TODO(Dominic): Consider joins.
   * TODO(Dominic): Potentially, there is still a bug if two neighbouring ranks are
   *                forking at the same time.
   *
   */
  bool hasToReceiveMetadata(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int fromRank);

  /**
   * Similar to hasToReceiveMetadata. However ignores
   * that the src rank in the adjacency information
   * might be a forking/joining one.
   */
  bool hasToReceiveMetadataIgnoreForksJoins(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int fromRank);



  /**
   * Checks for all cell descriptions (ADER-DG, FV, ...)
   * corresponding to the heap index at position \p src
   * in getCellDescriptions() if now is the time to
   * send out face data to a neighbouring rank
   *
   * The face corresponding to the adjacent cells \p src and \p dest
   * must be an inside face if periodic boundary conditions
   * are switched off.
   *
   * \note This method should only be used
   * if hasToSendMetadata(...) returns true.
   *
   * <h2>Periodic boundary conditions<h2>
   * If periodic boundary conditions are switched on,
   * fhe face corresponding to the adjacent cells
   * \p src and \src dest might be inside, outside,
   * or on the boundary.
   *
   * <h2>Face data exchange counters<\h2>
   * On every cell description, we hold a field of 2*d
   * counters. If a face is part of the MPI boundary,
   * we initialise the corresponding counter with
   * value 2^{d-1}.
   *
   * In the Prediction::prepareSendToNeighbour(...) and
   * RiemannSolver::mergeWithNeighbour(...) routine,
   * we then decrement the counters for the face
   * every time one of the 2^{d-1}
   * adjacent vertices touches the face.
   *
   * \see decrementCounters
   */
  bool hasToSendDataToNeighbour(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest) const;

  /**
   * Checks for all cell descriptions (ADER-DG, FV, ...)
   * corresponding to the heap index at position \p dest
   * in getCellDescriptions() if now is the time to
   * receive face data from a neighbouring rank
   *
   * The face corresponding to the adjacent cells \p dest and \p src
   * must be an inside face if periodic boundary conditions
   * are switched off.
   *
   * \note This method should only be used
   * if hasToReceiveMetadata(...) returns true.
   *
   * <h2>Periodic boundary conditions<h2>
   * If periodic boundary conditions are switched on,
   * fhe face corresponding to the adjacent cells
   * \p src and \src dest might be inside, outside,
   * or on the boundary.
   *
   * <h2>Face data exchange counters<\h2>
   * On every cell description, we hold a field of 2*d
   * counters. If a face is part of the MPI boundary,
   * we initialise the corresponding counter with
   * value 2^{d-1}.
   *
   * In the Prediction::prepareSendToNeighbour(...) and
   * RiemannSolver::mergeWithNeighbour(...) routine,
   * we then decrement the counters for the face
   * every time one of the 2^{d-1}
   * adjacent vertices touches the face.
   *
   * \see decrementCounters
   */
  bool hasToMergeWithNeighbourData(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest) const;

  /**
   * Every call of this function decrements the
   * faceDataExchangeCounter for the face corresponding
   * to the source and destination position pair \p src and \p dest
   * for all cell descriptions corresponding to \p cellDescriptionsIndex.
   *
   * \note Unfortunately, we cannot move the face data exchange
   * counters from the cell descriptions (ADERDGCellDescription,
   * FiniteVolumesCellDescription,...)
   * to the fineGridCell since we access the cell descriptions
   * from the vertices in Prediction::prepareToSendToNeighbour(...)
   * and RiemannSolver::mergeWithNeighbour(...). We do not have access
   * to the corresponding cell records in this case.
   *
   * Naturally, we cannot use counters if we do not
   * have any cell description registered on a cell
   * adjacent to a vertex. In this case,
   * the function simply returns.
   *
   * \see hasToSendFace
   */
  void tryDecrementFaceDataExchangeCountersOfSource(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest) const;

  /**
   * Resets the face data exchange coutnesr of
   * the the cell descriptions corresponding
   * to cell position \p dest.
   *
   * \note Requires that hasToReceiveSolverData(...) has returned true
   * beforehand.
   *
   * \see hasToMergeNeighbourData
   */
  void setFaceDataExchangeCountersOfDestination(
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const int value) const;

#endif
};

#endif
