#include "exahype/records/RepositoryState.h"

exahype::records::RepositoryState::PersistentRecords::PersistentRecords() {
   
}


exahype::records::RepositoryState::PersistentRecords::PersistentRecords(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices):
_action(action),
_numberOfIterations(numberOfIterations),
_exchangeBoundaryVertices(exchangeBoundaryVertices) {
   
}

exahype::records::RepositoryState::RepositoryState() {
   
}


exahype::records::RepositoryState::RepositoryState(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._action, persistentRecords._numberOfIterations, persistentRecords._exchangeBoundaryVertices) {
   
}


exahype::records::RepositoryState::RepositoryState(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices):
_persistentRecords(action, numberOfIterations, exchangeBoundaryVertices) {
   
}


exahype::records::RepositoryState::~RepositoryState() { }

std::string exahype::records::RepositoryState::toString(const Action& param) {
   switch (param) {
      case WriteCheckpoint: return "WriteCheckpoint";
      case ReadCheckpoint: return "ReadCheckpoint";
      case Terminate: return "Terminate";
      case RunOnAllNodes: return "RunOnAllNodes";
      case UseAdapterMeshRefinement: return "UseAdapterMeshRefinement";
      case UseAdapterPlotAugmentedAMRGrid: return "UseAdapterPlotAugmentedAMRGrid";
      case UseAdapterInitialConditionAndTimeStepSizeComputation: return "UseAdapterInitialConditionAndTimeStepSizeComputation";
      case UseAdapterPredictionAndFusedTimeSteppingInitialisation: return "UseAdapterPredictionAndFusedTimeSteppingInitialisation";
      case UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot: return "UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot";
      case UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d: return "UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d";
      case UseAdapterGridErasing: return "UseAdapterGridErasing";
      case UseAdapterADERDGTimeStep: return "UseAdapterADERDGTimeStep";
      case UseAdapterPlotAndADERDGTimeStep: return "UseAdapterPlotAndADERDGTimeStep";
      case UseAdapterPredictionRerun: return "UseAdapterPredictionRerun";
      case UseAdapterLimiterStatusSpreading: return "UseAdapterLimiterStatusSpreading";
      case UseAdapterLimiterStatusMergingAndSpreadingMPI: return "UseAdapterLimiterStatusMergingAndSpreadingMPI";
      case UseAdapterLimiterStatusMergingMPI: return "UseAdapterLimiterStatusMergingMPI";
      case UseAdapterReinitialisation: return "UseAdapterReinitialisation";
      case UseAdapterSolutionRecomputationAndTimeStepSizeComputation: return "UseAdapterSolutionRecomputationAndTimeStepSizeComputation";
      case UseAdapterNeighbourDataMerging: return "UseAdapterNeighbourDataMerging";
      case UseAdapterSolutionUpdate: return "UseAdapterSolutionUpdate";
      case UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation: return "UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation";
      case UseAdapterTimeStepSizeComputation: return "UseAdapterTimeStepSizeComputation";
      case UseAdapterPrediction: return "UseAdapterPrediction";
      case UseAdapterPredictionAndPlot: return "UseAdapterPredictionAndPlot";
      case UseAdapterPredictionAndPlot2d: return "UseAdapterPredictionAndPlot2d";
      case NumberOfAdapters: return "NumberOfAdapters";
   }
   return "undefined";
}

std::string exahype::records::RepositoryState::getActionMapping() {
   return "Action(WriteCheckpoint=0,ReadCheckpoint=1,Terminate=2,RunOnAllNodes=3,UseAdapterMeshRefinement=4,UseAdapterPlotAugmentedAMRGrid=5,UseAdapterInitialConditionAndTimeStepSizeComputation=6,UseAdapterPredictionAndFusedTimeSteppingInitialisation=7,UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot=8,UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d=9,UseAdapterGridErasing=10,UseAdapterADERDGTimeStep=11,UseAdapterPlotAndADERDGTimeStep=12,UseAdapterPredictionRerun=13,UseAdapterLimiterStatusSpreading=14,UseAdapterLimiterStatusMergingAndSpreadingMPI=15,UseAdapterLimiterStatusMergingMPI=16,UseAdapterReinitialisation=17,UseAdapterSolutionRecomputationAndTimeStepSizeComputation=18,UseAdapterNeighbourDataMerging=19,UseAdapterSolutionUpdate=20,UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation=21,UseAdapterTimeStepSizeComputation=22,UseAdapterPrediction=23,UseAdapterPredictionAndPlot=24,UseAdapterPredictionAndPlot2d=25,NumberOfAdapters=26)";
}


std::string exahype::records::RepositoryState::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::RepositoryState::toString (std::ostream& out) const {
   out << "("; 
   out << "action:" << toString(getAction());
   out << ",";
   out << "numberOfIterations:" << getNumberOfIterations();
   out << ",";
   out << "exchangeBoundaryVertices:" << getExchangeBoundaryVertices();
   out <<  ")";
}


exahype::records::RepositoryState::PersistentRecords exahype::records::RepositoryState::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::RepositoryStatePacked exahype::records::RepositoryState::convert() const{
   return RepositoryStatePacked(
      getAction(),
      getNumberOfIterations(),
      getExchangeBoundaryVertices()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::RepositoryState::_log( "exahype::records::RepositoryState" );
   
   MPI_Datatype exahype::records::RepositoryState::Datatype = 0;
   MPI_Datatype exahype::records::RepositoryState::FullDatatype = 0;
   
   
   void exahype::records::RepositoryState::initDatatype() {
      {
         RepositoryState dummyRepositoryState[2];
         
         const int Attributes = 4;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //action
            MPI_INT,		 //numberOfIterations
            MPI_CHAR,		 //exchangeBoundaryVertices
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //action
            1,		 //numberOfIterations
            1,		 //exchangeBoundaryVertices
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._action))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._numberOfIterations))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._exchangeBoundaryVertices))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[1]._persistentRecords._action))), 		&disp[3] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RepositoryState::Datatype );
         MPI_Type_commit( &RepositoryState::Datatype );
         
      }
      {
         RepositoryState dummyRepositoryState[2];
         
         const int Attributes = 4;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //action
            MPI_INT,		 //numberOfIterations
            MPI_CHAR,		 //exchangeBoundaryVertices
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //action
            1,		 //numberOfIterations
            1,		 //exchangeBoundaryVertices
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._action))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._numberOfIterations))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[0]._persistentRecords._exchangeBoundaryVertices))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryState[1]._persistentRecords._action))), 		&disp[3] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RepositoryState::FullDatatype );
         MPI_Type_commit( &RepositoryState::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::RepositoryState::shutdownDatatype() {
      MPI_Type_free( &RepositoryState::Datatype );
      MPI_Type_free( &RepositoryState::FullDatatype );
      
   }
   
   void exahype::records::RepositoryState::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      _senderDestinationRank = destination;
      
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::RepositoryState "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         
      }
      else {
      
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Isend(
            this, 1, Datatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      else {
         result = MPI_Isend(
            this, 1, FullDatatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message exahype::records::RepositoryState "
         << toString()
         << " to node " << destination
         << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "send(int)",msg.str() );
      }
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished send task for exahype::records::RepositoryState "
            << toString()
            << " sent to node " << destination
            << " failed: " << tarch::parallel::MPIReturnValueToString(result);
            _log.error("send(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::RepositoryState",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::RepositoryState",
            "send(int)", destination,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
         usleep(communicateSleep);
         
      }
      
      delete sendRequestHandle;
      #ifdef Debug
      _log.debug("send(int,int)", "sent " + toString() );
      #endif
      
   }
   
}



void exahype::records::RepositoryState::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   if (communicateSleep<0) {
   
      MPI_Status  status;
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
      _senderDestinationRank = status.MPI_SOURCE;
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive exahype::records::RepositoryState from node "
         << source << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "receive(int)", msg.str() );
      }
      
   }
   else {
   
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Irecv(
            this, 1, Datatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      else {
         result = MPI_Irecv(
            this, 1, FullDatatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive exahype::records::RepositoryState from node "
         << source << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "receive(int)", msg.str() );
      }
      
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished receive task for exahype::records::RepositoryState failed: "
            << tarch::parallel::MPIReturnValueToString(result);
            _log.error("receive(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "exahype::records::RepositoryState",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::RepositoryState",
            "receive(int)", source,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
         usleep(communicateSleep);
         
      }
      
      delete sendRequestHandle;
      
      _senderDestinationRank = status.MPI_SOURCE;
      #ifdef Debug
      _log.debug("receive(int,int)", "received " + toString() ); 
      #endif
      
   }
   
}



bool exahype::records::RepositoryState::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
   MPI_Status status;
   int  flag        = 0;
   MPI_Iprobe(
      MPI_ANY_SOURCE, tag,
      tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
   );
   if (flag) {
      int  messageCounter;
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         MPI_Get_count(&status, Datatype, &messageCounter);
      }
      else {
         MPI_Get_count(&status, FullDatatype, &messageCounter);
      }
      return messageCounter > 0;
   }
   else return false;
   
}

int exahype::records::RepositoryState::getSenderRank() const {
   assertion( _senderDestinationRank!=-1 );
   return _senderDestinationRank;
   
}
#endif


exahype::records::RepositoryStatePacked::PersistentRecords::PersistentRecords() {

}


exahype::records::RepositoryStatePacked::PersistentRecords::PersistentRecords(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices):
_action(action),
_numberOfIterations(numberOfIterations),
_exchangeBoundaryVertices(exchangeBoundaryVertices) {

}

exahype::records::RepositoryStatePacked::RepositoryStatePacked() {

}


exahype::records::RepositoryStatePacked::RepositoryStatePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._action, persistentRecords._numberOfIterations, persistentRecords._exchangeBoundaryVertices) {

}


exahype::records::RepositoryStatePacked::RepositoryStatePacked(const Action& action, const int& numberOfIterations, const bool& exchangeBoundaryVertices):
_persistentRecords(action, numberOfIterations, exchangeBoundaryVertices) {

}


exahype::records::RepositoryStatePacked::~RepositoryStatePacked() { }

std::string exahype::records::RepositoryStatePacked::toString(const Action& param) {
return exahype::records::RepositoryState::toString(param);
}

std::string exahype::records::RepositoryStatePacked::getActionMapping() {
return exahype::records::RepositoryState::getActionMapping();
}



std::string exahype::records::RepositoryStatePacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::RepositoryStatePacked::toString (std::ostream& out) const {
out << "("; 
out << "action:" << toString(getAction());
out << ",";
out << "numberOfIterations:" << getNumberOfIterations();
out << ",";
out << "exchangeBoundaryVertices:" << getExchangeBoundaryVertices();
out <<  ")";
}


exahype::records::RepositoryStatePacked::PersistentRecords exahype::records::RepositoryStatePacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::RepositoryState exahype::records::RepositoryStatePacked::convert() const{
return RepositoryState(
   getAction(),
   getNumberOfIterations(),
   getExchangeBoundaryVertices()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::RepositoryStatePacked::_log( "exahype::records::RepositoryStatePacked" );

MPI_Datatype exahype::records::RepositoryStatePacked::Datatype = 0;
MPI_Datatype exahype::records::RepositoryStatePacked::FullDatatype = 0;


void exahype::records::RepositoryStatePacked::initDatatype() {
   {
      RepositoryStatePacked dummyRepositoryStatePacked[2];
      
      const int Attributes = 4;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //action
         MPI_INT,		 //numberOfIterations
         MPI_CHAR,		 //exchangeBoundaryVertices
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //action
         1,		 //numberOfIterations
         1,		 //exchangeBoundaryVertices
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._action))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._numberOfIterations))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._exchangeBoundaryVertices))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[1]._persistentRecords._action))), 		&disp[3] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RepositoryStatePacked::Datatype );
      MPI_Type_commit( &RepositoryStatePacked::Datatype );
      
   }
   {
      RepositoryStatePacked dummyRepositoryStatePacked[2];
      
      const int Attributes = 4;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //action
         MPI_INT,		 //numberOfIterations
         MPI_CHAR,		 //exchangeBoundaryVertices
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //action
         1,		 //numberOfIterations
         1,		 //exchangeBoundaryVertices
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._action))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._numberOfIterations))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[0]._persistentRecords._exchangeBoundaryVertices))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRepositoryStatePacked[1]._persistentRecords._action))), 		&disp[3] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RepositoryStatePacked::FullDatatype );
      MPI_Type_commit( &RepositoryStatePacked::FullDatatype );
      
   }
   
}


void exahype::records::RepositoryStatePacked::shutdownDatatype() {
   MPI_Type_free( &RepositoryStatePacked::Datatype );
   MPI_Type_free( &RepositoryStatePacked::FullDatatype );
   
}

void exahype::records::RepositoryStatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   _senderDestinationRank = destination;
   
   if (communicateSleep<0) {
   
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message exahype::records::RepositoryStatePacked "
         << toString()
         << " to node " << destination
         << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "send(int)",msg.str() );
      }
      
   }
   else {
   
   MPI_Request* sendRequestHandle = new MPI_Request();
   MPI_Status   status;
   int          flag = 0;
   int          result;
   
   clock_t      timeOutWarning   = -1;
   clock_t      timeOutShutdown  = -1;
   bool         triggeredTimeoutWarning = false;
   
   if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Isend(
         this, 1, Datatype, destination,
         tag, tarch::parallel::Node::getInstance().getCommunicator(),
         sendRequestHandle
      );
      
   }
   else {
      result = MPI_Isend(
         this, 1, FullDatatype, destination,
         tag, tarch::parallel::Node::getInstance().getCommunicator(),
         sendRequestHandle
      );
      
   }
   if  (result!=MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message exahype::records::RepositoryStatePacked "
      << toString()
      << " to node " << destination
      << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error( "send(int)",msg.str() );
   }
   result = MPI_Test( sendRequestHandle, &flag, &status );
   while (!flag) {
      if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test( sendRequestHandle, &flag, &status );
      if (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "testing for finished send task for exahype::records::RepositoryStatePacked "
         << toString()
         << " sent to node " << destination
         << " failed: " << tarch::parallel::MPIReturnValueToString(result);
         _log.error("send(int)", msg.str() );
      }
      
      // deadlock aspect
      if (
         tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
         (clock()>timeOutWarning) &&
         (!triggeredTimeoutWarning)
      ) {
         tarch::parallel::Node::getInstance().writeTimeOutWarning(
         "exahype::records::RepositoryStatePacked",
         "send(int)", destination,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::RepositoryStatePacked",
         "send(int)", destination,tag,1
         );
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
      
   }
   
   delete sendRequestHandle;
   #ifdef Debug
   _log.debug("send(int,int)", "sent " + toString() );
   #endif
   
}

}



void exahype::records::RepositoryStatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

   MPI_Status  status;
   const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
   _senderDestinationRank = status.MPI_SOURCE;
   if ( result != MPI_SUCCESS ) {
      std::ostringstream msg;
      msg << "failed to start to receive exahype::records::RepositoryStatePacked from node "
      << source << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error( "receive(int)", msg.str() );
   }
   
}
else {

   MPI_Request* sendRequestHandle = new MPI_Request();
   MPI_Status   status;
   int          flag = 0;
   int          result;
   
   clock_t      timeOutWarning   = -1;
   clock_t      timeOutShutdown  = -1;
   bool         triggeredTimeoutWarning = false;
   
   if (exchangeOnlyAttributesMarkedWithParallelise) {
      result = MPI_Irecv(
         this, 1, Datatype, source, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
      );
      
   }
   else {
      result = MPI_Irecv(
         this, 1, FullDatatype, source, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
      );
      
   }
   if ( result != MPI_SUCCESS ) {
      std::ostringstream msg;
      msg << "failed to start to receive exahype::records::RepositoryStatePacked from node "
      << source << ": " << tarch::parallel::MPIReturnValueToString(result);
      _log.error( "receive(int)", msg.str() );
   }
   
   result = MPI_Test( sendRequestHandle, &flag, &status );
   while (!flag) {
      if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
      if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
      result = MPI_Test( sendRequestHandle, &flag, &status );
      if (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "testing for finished receive task for exahype::records::RepositoryStatePacked failed: "
         << tarch::parallel::MPIReturnValueToString(result);
         _log.error("receive(int)", msg.str() );
      }
      
      // deadlock aspect
      if (
         tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
         (clock()>timeOutWarning) &&
         (!triggeredTimeoutWarning)
      ) {
         tarch::parallel::Node::getInstance().writeTimeOutWarning(
         "exahype::records::RepositoryStatePacked",
         "receive(int)", source,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::RepositoryStatePacked",
         "receive(int)", source,tag,1
         );
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
      usleep(communicateSleep);
      
   }
   
   delete sendRequestHandle;
   
   _senderDestinationRank = status.MPI_SOURCE;
   #ifdef Debug
   _log.debug("receive(int,int)", "received " + toString() ); 
   #endif
   
}

}



bool exahype::records::RepositoryStatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
MPI_Status status;
int  flag        = 0;
MPI_Iprobe(
   MPI_ANY_SOURCE, tag,
   tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
);
if (flag) {
   int  messageCounter;
   if (exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Get_count(&status, Datatype, &messageCounter);
   }
   else {
      MPI_Get_count(&status, FullDatatype, &messageCounter);
   }
   return messageCounter > 0;
}
else return false;

}

int exahype::records::RepositoryStatePacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif



