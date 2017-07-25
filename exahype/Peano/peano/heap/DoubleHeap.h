// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _PEANO_HEAP_DOUBLE_HEAP_H_
#define _PEANO_HEAP_DOUBLE_HEAP_H_


#include "peano/heap/Heap.h"


namespace peano {
  namespace heap {
    template<
      class MasterWorkerExchanger,
      class JoinForkExchanger,
      class NeighbourDataExchanger,
      // @tood Perhaps remove default and always align?
      class VectorContainer = std::vector<double>
    >
    class DoubleHeap;


    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      PlainBoundaryDataExchanger< double, true >
    >     PlainDoubleHeap;

    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      PlainBoundaryDataExchanger< double, true >,
      std::vector< double, HeapAllocator<double, 32 > >
    >     PlainDoubleHeapAlignment32;

    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      PlainBoundaryDataExchanger< double, true >,
      std::vector< double, HeapAllocator<double, 64 > >
    >     PlainDoubleHeapAlignment64;

    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      RLEBoundaryDataExchanger< double, true >
    >     RLEDoubleHeap;

    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      RLEBoundaryDataExchanger< double, true >,
      std::vector< double, HeapAllocator<double, 32 > >
    >     RLEDoubleHeapAlignment32;

    typedef DoubleHeap<
      SynchronousDataExchanger< double, true >,
      SynchronousDataExchanger< double, true >,
      RLEBoundaryDataExchanger< double, true >,
      std::vector< double, HeapAllocator<double, 64 > >
    >     RLEDoubleHeapAlignment64;

  }
}



/**
 * DoubleHeap
 *
 * This is a specialised variant of the heap for doubles. It works directly
 * with doubles held in a std::vector. It does not rely on DaStGen for the
 * data at all and not wrap any data into DaStGen records. It thus should be
 * faster than the standard version.
 *
 * <h1> Working with plain double pointer </h1>
 *
 * With this class, you may use getData().data() yielding a plain double
 * pointer. It is probably aligned if you choose alignment.
 *
 *
 * <h1> Alignment </h1>
 *
 * A big difference to the standard heap class is that this class can work with
 * aligned data structuures. This makes the class however incompatible with
 * other std::vector<double> instances where no alignment is used. Please consult
 * the HeapAllocator for defails on the alignment.
 *
 * <h1> Method documentation </h1>
 *
 * This is a specialisation of the general-purpose heap. As such, the
 * documentation for the routines here is empty. Please study peano::heap::Heap
 * to find out about the intended semantics.
 *
 *
 * @author Tobias Weinzierl
 */
template <class MasterWorkerExchanger, class JoinForkExchanger, class NeighbourDataExchanger, class VectorContainer>
class peano::heap::DoubleHeap: public tarch::services::Service, peano::heap::AbstractHeap {
  private:
    static tarch::logging::Log _log;

    typedef std::map<int, VectorContainer*>  HeapContainer;

    HeapContainer    _heapData;

    std::list<int>   _deletedHeapIndices;

    std::list<int>   _recycledHeapIndices;

    int _nextIndex;

    #ifdef Parallel
    int                                    _neighbourDataExchangerMetaDataTag;
    int                                    _neighbourDataExchangerDataTag;

    MasterWorkerExchanger                  _masterWorkerExchanger;
    JoinForkExchanger                      _joinForkExchanger;
    std::map<int, NeighbourDataExchanger>  _neighbourDataExchanger;
    #endif

    int _maximumNumberOfHeapEntries;

    int _numberOfHeapAllocations;

    int _numberOfHeapFrees;

    std::string _name;

    DoubleHeap();

    ~DoubleHeap();

  public:
    enum class Allocation {
      DoNotUseAnyRecycledEntry,
      UseOnlyRecycledEntries,
      UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired
    };

    typedef VectorContainer  HeapEntries;

    virtual void startToSendSynchronousData();

    virtual void startToSendBoundaryData(bool isTraversalInverted);

    virtual void finishedToSendSynchronousData();

    virtual void finishedToSendBoundaryData(bool isTraversalInverted);

    static DoubleHeap& getInstance();

    HeapEntries& getData(int index);

    const HeapEntries& getData(int index) const;

    int createData(int numberOfEntries=0, int initialCapacity=0, Allocation allocation = Allocation::UseRecycledEntriesIfPossibleCreateNewEntriesIfRequired);

    void createDataForIndex(int wantedIndex, int numberOfEntries=0, int initialCapacity=0);

    void reserveHeapEntriesForRecycling(int numberOfEntries);

    bool areRecycleEntriesAvailable() const;

    bool isValidIndex(int index) const;

    void deleteData(int index, bool recycle = false);

    void deleteAllData();

    int getNumberOfAllocatedEntries() const;

    void moveData( int toIndex, int fromIndex );

    void addData( int index, const HeapEntries& entries );

    void addData( int index, const double&        entry );

    void restart();

    void shutdown();

    void setName(std::string name);

    void sendData(
      int                                           index,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    void sendData(
      const HeapEntries&                            data,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    void sendData(
      const double*                                 data,
      int                                           size,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    HeapEntries receiveData(
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    void receiveData(
      double*                                       data,
      int                                           size,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    int receiveData(
      int                                           index,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  position,
      int                                           level,
      MessageType                                   messageType
    );

    virtual void receiveDanglingMessages();

    std::string toString() const;

    void plotStatistics() const;

    void clearStatistics();

    void logContentToWarningDevice();
};



#include "peano/heap/DoubleHeap.cpph"


#endif
