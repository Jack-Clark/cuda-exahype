#include "mpibalancing/FairNodePoolStrategy.h"
#include "tarch/parallel/Node.h"
#include "tarch/compiler/CompilerSpecificSettings.h"


#include <cstdlib>

#include <sstream>
#include <limits>


tarch::logging::Log mpibalancing::FairNodePoolStrategy::_log( "mpibalancing::FairNodePoolStrategy" );


mpibalancing::FairNodePoolStrategy::FairNodePoolStrategy(int mpiRanksPerNode, double waitTimeOutSec):
  NodePoolStrategy(),
  _tag(-1),
  _nodes(),
  _waitTimeOut(waitTimeOutSec),
  _mpiRanksPerNode(mpiRanksPerNode) {

  if (
    tarch::parallel::Node::getInstance().getNumberOfNodes() <= _mpiRanksPerNode &&
    tarch::parallel::Node::getInstance().isGlobalMaster()
  ) {
    logWarning(
      "FairNodePoolStrategy(int,double)",
      "number of mpi ranks is smaller than modulo argument of FairNodePoolStrategy. First " <<
      _mpiRanksPerNode << " will not be assigned any job, so setting reduces to serial run"
    );
  }
}


mpibalancing::FairNodePoolStrategy::~FairNodePoolStrategy() {
}


void mpibalancing::FairNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  #ifdef Parallel
  assertion( _tag >= 0 );

  const std::clock_t waitTimeoutTimeStamp = clock() + static_cast<std::clock_t>(std::floor(_waitTimeOut * CLOCKS_PER_SEC));

  bool continueToWait = true;
  while ( continueToWait ) {
    while (tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_tag, true)) {
      tarch::parallel::messages::WorkerRequestMessage message;
      message.receive(MPI_ANY_SOURCE,_tag, true, SendAndReceiveLoadBalancingMessagesBlocking);
      queue.push_back( message );
    }

    continueToWait =
        hasIdleNode(-1) &&
        (static_cast<int>(queue.size())+1 < getNumberOfRegisteredNodes()-getNumberOfIdleNodes()) &&
        (clock() < waitTimeoutTimeStamp);
  }
  #endif
}


void mpibalancing::FairNodePoolStrategy::logQueue( const RequestQueue& queue ) const {
  #ifdef Parallel
  if (queue.empty()) {
  _log.debug( "logQueue()", "queue is empty" );
  }
  else {
    std::ostringstream msg;
    msg << "queue: ";

    for (RequestQueue::const_iterator p = queue.begin(); p != queue.end(); p++ ) {
      msg << p->getSenderRank() << ",";
    }
    _log.debug( "logQueue()", msg.str() );
  }
  #endif
}


tarch::parallel::messages::WorkerRequestMessage mpibalancing::FairNodePoolStrategy::extractElementFromRequestQueue(RequestQueue& queue) {
  assertion( !queue.empty() );

  RequestQueue::iterator pResultInQueue;

  #ifdef Parallel
  int                    workersOfRankCurrentlyAnswered = std::numeric_limits<int>::max();
  for (RequestQueue::iterator p = queue.begin(); p != queue.end(); p++) {
    if (getWorkersOfNode(p->getSenderRank()) < workersOfRankCurrentlyAnswered) {
      workersOfRankCurrentlyAnswered = getWorkersOfNode(p->getSenderRank());
      pResultInQueue                 = p;
    }
  }
  #endif

  tarch::parallel::messages::WorkerRequestMessage result = *pResultInQueue;
  queue.erase(pResultInQueue);

  return result;
}


int mpibalancing::FairNodePoolStrategy::getWorkersOfNode( int rank ) const {
  int result = -1;

  for (NodeContainer::const_iterator p=_nodes.begin(); p!=_nodes.end(); p++) {
    if (p->getRank()==rank) {
      result = p->getNumberOfBookedWorkers();
    }
  }

  assertion( result>=0 );
  return result;
}


void mpibalancing::FairNodePoolStrategy::addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node) {
  #ifdef Parallel
  assertion( !isRegisteredNode(node.getSenderRank()) );

  logTraceInWith1Argument( "addNode(...)", node.getSenderRank() );
  NodePoolListEntry newEntry(
    _mpiRanksPerNode,
    node.getSenderRank(),
    tarch::parallel::StringTools::convert(node.getNodeName())
  );
  _nodes.push_back( newEntry ) ;
  _nodes.sort();
  logTraceOutWith1Argument( "addNode(...)", newEntry.toString() );
  #endif
}


void mpibalancing::FairNodePoolStrategy::removeNode( int rank ) {
  assertion( isRegisteredNode(rank) );

  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      #ifdef Debug
      _log.debug( "removeNode(int)", "remove entry " + p->toString() );
      #endif
      _nodes.erase(p);
      _nodes.sort();
      return;
    }
  }
}


bool mpibalancing::FairNodePoolStrategy::hasIdleNode(int forMaster) const {
  if (_nodes.empty()) return false;

  if (forMaster==AnyMaster) {
    return _nodes.front().isIdle() || _nodes.back().isIdle();
  }
  else {
    return _nodes.front().isIdle();
  }
}


int mpibalancing::FairNodePoolStrategy::removeNextIdleNode() {
  assertion1( !_nodes.empty(), _nodes.size() );
  assertion4(
    FairNodePoolStrategy::hasIdleNode(NodePoolStrategy::AnyMaster), _nodes.size(),
    _nodes.front().toString(), _nodes.back().toString(), toString()
  );


  int result = -1;

  if (_nodes.back().isIdle()) {
    result = _nodes.back().getRank();
    _nodes.pop_back();
  }
  else
  if (_nodes.front().isIdle()) {
    result = _nodes.front().getRank();
    _nodes.pop_front();
  }

  _nodes.sort();

  return result;
}


int mpibalancing::FairNodePoolStrategy::getNumberOfIdleNodes() const {
  int result = 0;
  NodeContainer::const_iterator p = _nodes.begin();
  while (p != _nodes.end()&& p->isIdle() ) {
    p++;
    result++;
  }
  return result;
}


void mpibalancing::FairNodePoolStrategy::setNodeIdle( int rank ) {
  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->deActivate();
    }
  }

  _nodes.sort();
}


bool mpibalancing::FairNodePoolStrategy::isRegisteredNode(int rank) const {
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      return true;
    }
  }
  return false;
}


bool mpibalancing::FairNodePoolStrategy::isIdleNode(int rank) const {
  assertion1( isRegisteredNode(rank), rank );
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank && p->isIdle() ) {
      return true;
    }
  }
  return false;
}


int mpibalancing::FairNodePoolStrategy::getNumberOfRegisteredNodes() const {
  return static_cast<int>( _nodes.size() );
}


std::string mpibalancing::FairNodePoolStrategy::toString() const {
  std::ostringstream result;
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    p->toString(result);
  }
  return result.str();
}


void mpibalancing::FairNodePoolStrategy::reserveParticularNode(int rank) {
  assertion( rank>=0 );
  assertion1( isIdleNode(rank), rank );

  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->activate();
    }
  }

  _nodes.sort();
}


int mpibalancing::FairNodePoolStrategy::reserveNode(int forMaster) {
  assertion1(hasIdleNode(forMaster),forMaster);

  NodePoolListEntry result;

  if (_nodes.front().isIdle()) {
    result = _nodes.front();
    _nodes.pop_front();
  }
  else {
    result = _nodes.back();
    _nodes.pop_back();
  }

  result.activate();
  _nodes.push_back(result);
  _nodes.sort();

  for (NodeContainer::iterator p=_nodes.begin(); p!=_nodes.end(); p++ ) {
    if (p->getRank()==forMaster) {
      p->addNewWorker();
    }
  }

  updateNodeWeights();

  return result.getRank();
}


void mpibalancing::FairNodePoolStrategy::updateNodeWeights() {
  for (NodeContainer::iterator p=_nodes.begin(); p!=_nodes.end(); p++ ) {
    p->reduceNumberOfBookedWorkers();
  }
}


void mpibalancing::FairNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}



mpibalancing::FairNodePoolStrategy::NodePoolListEntry::NodePoolListEntry():
  _mpiRanksPerNode(-1),
  _rank(-1),
  _bookedWorkers(-1),
  _state(WORKING),
  _name("undef") {
}


mpibalancing::FairNodePoolStrategy::NodePoolListEntry::NodePoolListEntry(
  const int   mpiRanksPerNode,
  int         rank,
  const std::string& name
  ):
  _mpiRanksPerNode(mpiRanksPerNode),
  _rank(rank),
  _bookedWorkers(0),
  _state(WORKING),
  _name(name) {
}


mpibalancing::FairNodePoolStrategy::NodePoolListEntry::~NodePoolListEntry() {
}


std::string mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getNodeName() const {
  return _name;
}


std::string mpibalancing::FairNodePoolStrategy::NodePoolListEntry::toString() const {
  std::ostringstream out;
  toString(out);
  return out.str();
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::toString(std::ostream& out) const {
  out << "(rank:" << _rank;
  switch (_state) {
    case IDLE:    out << ",state:idle";     break;
    case WORKING: out << ",state:working";  break;
  }
  out << ",booked-workers:" << _bookedWorkers;
  out << ",name:" << _name << ")";
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::operator==( const NodePoolListEntry& than ) const {
  return _rank==than._rank;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::activate() {
  assertionEquals1( _state, IDLE, toString() );
  _state = WORKING;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::deActivate() {
  _state         = IDLE;
  _bookedWorkers = 0;
}


int mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getRank() const {
  return _rank;
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::isIdle() const {
  return _state == IDLE;
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::operator<( const mpibalancing::FairNodePoolStrategy::NodePoolListEntry& than ) const {
  const bool lhsIsClassA = isIdle()       && _rank>=_mpiRanksPerNode;
  const bool rhsIsClassA = than.isIdle()  && than._rank>=_mpiRanksPerNode;
  const bool lhsIsClassB = !isIdle()      && _rank>=_mpiRanksPerNode;
  const bool rhsIsClassB = !than.isIdle() && than._rank>=_mpiRanksPerNode;
  const bool lhsIsClassC = _rank<_mpiRanksPerNode;
  const bool rhsIsClassC = than._rank<_mpiRanksPerNode;

  const bool smallerDueToClassOrder =
    (_rank % _mpiRanksPerNode <  than._rank % _mpiRanksPerNode) ||
    (_rank % _mpiRanksPerNode == than._rank % _mpiRanksPerNode && _rank < than._rank);

  return
         (lhsIsClassA && !rhsIsClassA)
      || (lhsIsClassB &&  rhsIsClassC)
      || (lhsIsClassA &&  rhsIsClassA && smallerDueToClassOrder)
      || (lhsIsClassB &&  rhsIsClassB && smallerDueToClassOrder)
      || (lhsIsClassC &&  rhsIsClassC && smallerDueToClassOrder);
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::addNewWorker() {
  _bookedWorkers+=1.0;
}


int mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getNumberOfBookedWorkers() const {
  return _bookedWorkers;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::reduceNumberOfBookedWorkers() {
  _bookedWorkers /= 2.0;
  assertion( _bookedWorkers>=0.0 );
}
