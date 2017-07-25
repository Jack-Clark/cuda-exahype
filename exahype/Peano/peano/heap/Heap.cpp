#include "peano/heap/Heap.h"


std::string peano::heap::toString(MessageType type) {
  switch (type) {
    case MessageType::NeighbourCommunication:
      return "neighbour";
    case MessageType::ForkOrJoinCommunication:
      return "fork-or-join";
    case MessageType::MasterWorkerCommunication:
      return "master-worker";
  }

  assertionMsg( false, "should not be entered ever" );
  return "<undef>";
}
