#include "exahype/records/State.h"

#if !defined(TrackGridStatistics) && !defined(Parallel)
   exahype::records::State::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::State::PersistentRecords::PersistentRecords(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _mergeMode(mergeMode),
   _sendMode(sendMode),
   _reinitTimeStepData(reinitTimeStepData),
   _stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
   _timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
   _limiterDomainHasChanged(limiterDomainHasChanged),
   _hasRefined(hasRefined),
   _hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
   _hasErased(hasErased),
   _hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
   _hasChangedVertexOrCellState(hasChangedVertexOrCellState),
   _hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
   _isTraversalInverted(isTraversalInverted) {
      
   }
   
   exahype::records::State::State() {
      
   }
   
   
   exahype::records::State::State(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted) {
      
   }
   
   
   exahype::records::State::State(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _persistentRecords(mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
      
   }
   
   
   exahype::records::State::~State() { }
   
   std::string exahype::records::State::toString(const MergeMode& param) {
      switch (param) {
         case MergeNothing: return "MergeNothing";
         case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
         case MergeFaceData: return "MergeFaceData";
         case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
      }
      return "undefined";
   }
   
   std::string exahype::records::State::getMergeModeMapping() {
      return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,BroadcastAndMergeTimeStepDataAndMergeFaceData=3)";
   }
   std::string exahype::records::State::toString(const SendMode& param) {
      switch (param) {
         case SendNothing: return "SendNothing";
         case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
         case SendFaceData: return "SendFaceData";
         case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
      }
      return "undefined";
   }
   
   std::string exahype::records::State::getSendModeMapping() {
      return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
   }
   
   
   std::string exahype::records::State::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::State::toString (std::ostream& out) const {
      out << "("; 
      out << "mergeMode:" << toString(getMergeMode());
      out << ",";
      out << "sendMode:" << toString(getSendMode());
      out << ",";
      out << "reinitTimeStepData:" << getReinitTimeStepData();
      out << ",";
      out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
      out << ",";
      out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
      out << ",";
      out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
      out << ",";
      out << "hasRefined:" << getHasRefined();
      out << ",";
      out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
      out << ",";
      out << "hasErased:" << getHasErased();
      out << ",";
      out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
      out << ",";
      out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
      out << ",";
      out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
      out << ",";
      out << "isTraversalInverted:" << getIsTraversalInverted();
      out <<  ")";
   }
   
   
   exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::StatePacked exahype::records::State::convert() const{
      return StatePacked(
         getMergeMode(),
         getSendMode(),
         getReinitTimeStepData(),
         getStabilityConditionOfOneSolverWasViolated(),
         getTimeStepSizeWeightForPredictionRerun(),
         getLimiterDomainHasChanged(),
         getHasRefined(),
         getHasTriggeredRefinementForNextIteration(),
         getHasErased(),
         getHasTriggeredEraseForNextIteration(),
         getHasChangedVertexOrCellState(),
         getHasModifiedGridInPreviousIteration(),
         getIsTraversalInverted()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );
      
      MPI_Datatype exahype::records::State::Datatype = 0;
      MPI_Datatype exahype::records::State::FullDatatype = 0;
      
      
      void exahype::records::State::initDatatype() {
         {
            State dummyState[2];
            
            const int Attributes = 11;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //mergeMode
               MPI_INT,		 //sendMode
               MPI_CHAR,		 //limiterDomainHasChanged
               MPI_CHAR,		 //hasRefined
               MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
               MPI_CHAR,		 //hasErased
               MPI_CHAR,		 //hasTriggeredEraseForNextIteration
               MPI_CHAR,		 //hasChangedVertexOrCellState
               MPI_CHAR,		 //hasModifiedGridInPreviousIteration
               MPI_CHAR,		 //isTraversalInverted
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //mergeMode
               1,		 //sendMode
               1,		 //limiterDomainHasChanged
               1,		 //hasRefined
               1,		 //hasTriggeredRefinementForNextIteration
               1,		 //hasErased
               1,		 //hasTriggeredEraseForNextIteration
               1,		 //hasChangedVertexOrCellState
               1,		 //hasModifiedGridInPreviousIteration
               1,		 //isTraversalInverted
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._mergeMode))), 		&disp[10] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype );
            MPI_Type_commit( &State::Datatype );
            
         }
         {
            State dummyState[2];
            
            const int Attributes = 14;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //mergeMode
               MPI_INT,		 //sendMode
               MPI_CHAR,		 //reinitTimeStepData
               MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
               MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
               MPI_CHAR,		 //limiterDomainHasChanged
               MPI_CHAR,		 //hasRefined
               MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
               MPI_CHAR,		 //hasErased
               MPI_CHAR,		 //hasTriggeredEraseForNextIteration
               MPI_CHAR,		 //hasChangedVertexOrCellState
               MPI_CHAR,		 //hasModifiedGridInPreviousIteration
               MPI_CHAR,		 //isTraversalInverted
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //mergeMode
               1,		 //sendMode
               1,		 //reinitTimeStepData
               1,		 //stabilityConditionOfOneSolverWasViolated
               1,		 //timeStepSizeWeightForPredictionRerun
               1,		 //limiterDomainHasChanged
               1,		 //hasRefined
               1,		 //hasTriggeredRefinementForNextIteration
               1,		 //hasErased
               1,		 //hasTriggeredEraseForNextIteration
               1,		 //hasChangedVertexOrCellState
               1,		 //hasModifiedGridInPreviousIteration
               1,		 //isTraversalInverted
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reinitTimeStepData))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._mergeMode))), 		&disp[13] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype );
            MPI_Type_commit( &State::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::State::shutdownDatatype() {
         MPI_Type_free( &State::Datatype );
         MPI_Type_free( &State::FullDatatype );
         
      }
      
      void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         _senderDestinationRank = destination;
         
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::State "
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
            msg << "was not able to send message exahype::records::State "
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
               msg << "testing for finished send task for exahype::records::State "
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
               "exahype::records::State",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::State",
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
   
   
   
   void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         _senderDestinationRank = status.MPI_SOURCE;
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::State from node "
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
            msg << "failed to start to receive exahype::records::State from node "
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
               msg << "testing for finished receive task for exahype::records::State failed: "
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
               "exahype::records::State",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::State",
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
   
   
   
   bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int exahype::records::State::getSenderRank() const {
      assertion( _senderDestinationRank!=-1 );
      return _senderDestinationRank;
      
   }
#endif


exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
   if ((6 >= (8 * sizeof(short int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((6 < (8 * sizeof(short int))));
   
}


exahype::records::StatePacked::PersistentRecords::PersistentRecords(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_isTraversalInverted(isTraversalInverted) {
   setHasRefined(hasRefined);
   setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
   setHasErased(hasErased);
   setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
   setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
   setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
   if ((6 >= (8 * sizeof(short int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((6 < (8 * sizeof(short int))));
   
}

exahype::records::StatePacked::StatePacked() {
   if ((6 >= (8 * sizeof(short int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((6 < (8 * sizeof(short int))));
   
}


exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted) {
   if ((6 >= (8 * sizeof(short int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((6 < (8 * sizeof(short int))));
   
}


exahype::records::StatePacked::StatePacked(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_persistentRecords(mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
   if ((6 >= (8 * sizeof(short int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((6 < (8 * sizeof(short int))));
   
}


exahype::records::StatePacked::~StatePacked() { }

std::string exahype::records::StatePacked::toString(const MergeMode& param) {
   return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getMergeModeMapping() {
   return exahype::records::State::getMergeModeMapping();
}

std::string exahype::records::StatePacked::toString(const SendMode& param) {
   return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getSendModeMapping() {
   return exahype::records::State::getSendModeMapping();
}



std::string exahype::records::StatePacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::StatePacked::toString (std::ostream& out) const {
   out << "("; 
   out << "mergeMode:" << toString(getMergeMode());
   out << ",";
   out << "sendMode:" << toString(getSendMode());
   out << ",";
   out << "reinitTimeStepData:" << getReinitTimeStepData();
   out << ",";
   out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
   out << ",";
   out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
   out << ",";
   out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
   out << ",";
   out << "hasRefined:" << getHasRefined();
   out << ",";
   out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
   out << ",";
   out << "hasErased:" << getHasErased();
   out << ",";
   out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
   out << ",";
   out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
   out << ",";
   out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
   out << ",";
   out << "isTraversalInverted:" << getIsTraversalInverted();
   out <<  ")";
}


exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::State exahype::records::StatePacked::convert() const{
   return State(
      getMergeMode(),
      getSendMode(),
      getReinitTimeStepData(),
      getStabilityConditionOfOneSolverWasViolated(),
      getTimeStepSizeWeightForPredictionRerun(),
      getLimiterDomainHasChanged(),
      getHasRefined(),
      getHasTriggeredRefinementForNextIteration(),
      getHasErased(),
      getHasTriggeredEraseForNextIteration(),
      getHasChangedVertexOrCellState(),
      getHasModifiedGridInPreviousIteration(),
      getIsTraversalInverted()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );
   
   MPI_Datatype exahype::records::StatePacked::Datatype = 0;
   MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;
   
   
   void exahype::records::StatePacked::initDatatype() {
      {
         StatePacked dummyStatePacked[2];
         
         const int Attributes = 6;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //mergeMode
            MPI_INT,		 //sendMode
            MPI_CHAR,		 //limiterDomainHasChanged
            MPI_CHAR,		 //isTraversalInverted
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //mergeMode
            1,		 //sendMode
            1,		 //limiterDomainHasChanged
            1,		 //isTraversalInverted
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._mergeMode))), 		&disp[5] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype );
         MPI_Type_commit( &StatePacked::Datatype );
         
      }
      {
         StatePacked dummyStatePacked[2];
         
         const int Attributes = 9;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //mergeMode
            MPI_INT,		 //sendMode
            MPI_CHAR,		 //reinitTimeStepData
            MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
            MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
            MPI_CHAR,		 //limiterDomainHasChanged
            MPI_CHAR,		 //isTraversalInverted
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //mergeMode
            1,		 //sendMode
            1,		 //reinitTimeStepData
            1,		 //stabilityConditionOfOneSolverWasViolated
            1,		 //timeStepSizeWeightForPredictionRerun
            1,		 //limiterDomainHasChanged
            1,		 //isTraversalInverted
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._reinitTimeStepData))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._mergeMode))), 		&disp[8] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype );
         MPI_Type_commit( &StatePacked::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::StatePacked::shutdownDatatype() {
      MPI_Type_free( &StatePacked::Datatype );
      MPI_Type_free( &StatePacked::FullDatatype );
      
   }
   
   void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      _senderDestinationRank = destination;
      
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::StatePacked "
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
         msg << "was not able to send message exahype::records::StatePacked "
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
            msg << "testing for finished send task for exahype::records::StatePacked "
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
            "exahype::records::StatePacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::StatePacked",
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



void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   if (communicateSleep<0) {
   
      MPI_Status  status;
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
      _senderDestinationRank = status.MPI_SOURCE;
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive exahype::records::StatePacked from node "
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
         msg << "failed to start to receive exahype::records::StatePacked from node "
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
            msg << "testing for finished receive task for exahype::records::StatePacked failed: "
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
            "exahype::records::StatePacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::StatePacked",
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



bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::StatePacked::getSenderRank() const {
   assertion( _senderDestinationRank!=-1 );
   return _senderDestinationRank;
   
}
#endif



#elif defined(TrackGridStatistics) && defined(Parallel)
exahype::records::State::PersistentRecords::PersistentRecords() {

}


exahype::records::State::PersistentRecords::PersistentRecords(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_gridConstructionState(gridConstructionState),
_firstGridSetupIteration(firstGridSetupIteration),
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_minMeshWidth(minMeshWidth),
_maxMeshWidth(maxMeshWidth),
_numberOfInnerVertices(numberOfInnerVertices),
_numberOfBoundaryVertices(numberOfBoundaryVertices),
_numberOfOuterVertices(numberOfOuterVertices),
_numberOfInnerCells(numberOfInnerCells),
_numberOfOuterCells(numberOfOuterCells),
_numberOfInnerLeafVertices(numberOfInnerLeafVertices),
_numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
_numberOfOuterLeafVertices(numberOfOuterLeafVertices),
_numberOfInnerLeafCells(numberOfInnerLeafCells),
_numberOfOuterLeafCells(numberOfOuterLeafCells),
_maxLevel(maxLevel),
_hasRefined(hasRefined),
_hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
_hasErased(hasErased),
_hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
_hasChangedVertexOrCellState(hasChangedVertexOrCellState),
_hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
_isTraversalInverted(isTraversalInverted),
_reduceStateAndCell(reduceStateAndCell),
_couldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag),
_subWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork) {

}

exahype::records::State::State() {

}


exahype::records::State::State(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._gridConstructionState, persistentRecords._firstGridSetupIteration, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted, persistentRecords._reduceStateAndCell, persistentRecords._couldNotEraseDueToDecompositionFlag, persistentRecords._subWorkerIsInvolvedInJoinOrFork) {

}


exahype::records::State::State(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_persistentRecords(gridConstructionState, firstGridSetupIteration, mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {

}


exahype::records::State::~State() { }

std::string exahype::records::State::toString(const GridConstructionState& param) {
switch (param) {
   case Default: return "Default";
   case Veto: return "Veto";
   case Aggressive: return "Aggressive";
}
return "undefined";
}

std::string exahype::records::State::getGridConstructionStateMapping() {
return "GridConstructionState(Default=0,Veto=1,Aggressive=2)";
}
std::string exahype::records::State::toString(const MergeMode& param) {
switch (param) {
   case MergeNothing: return "MergeNothing";
   case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
   case MergeFaceData: return "MergeFaceData";
   case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
}
return "undefined";
}

std::string exahype::records::State::getMergeModeMapping() {
return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,BroadcastAndMergeTimeStepDataAndMergeFaceData=3)";
}
std::string exahype::records::State::toString(const SendMode& param) {
switch (param) {
   case SendNothing: return "SendNothing";
   case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
   case SendFaceData: return "SendFaceData";
   case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
}
return "undefined";
}

std::string exahype::records::State::getSendModeMapping() {
return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
}


std::string exahype::records::State::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::State::toString (std::ostream& out) const {
out << "("; 
out << "gridConstructionState:" << toString(getGridConstructionState());
out << ",";
out << "firstGridSetupIteration:" << getFirstGridSetupIteration();
out << ",";
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
out << ",";
out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
out << ",";
out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
out << ",";
out << "numberOfInnerCells:" << getNumberOfInnerCells();
out << ",";
out << "numberOfOuterCells:" << getNumberOfOuterCells();
out << ",";
out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
out << ",";
out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
out << ",";
out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
out << ",";
out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
out << ",";
out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
out << ",";
out << "maxLevel:" << getMaxLevel();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out << ",";
out << "reduceStateAndCell:" << getReduceStateAndCell();
out << ",";
out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
out << ",";
out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
out <<  ")";
}


exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::StatePacked exahype::records::State::convert() const{
return StatePacked(
   getGridConstructionState(),
   getFirstGridSetupIteration(),
   getMergeMode(),
   getSendMode(),
   getReinitTimeStepData(),
   getStabilityConditionOfOneSolverWasViolated(),
   getTimeStepSizeWeightForPredictionRerun(),
   getLimiterDomainHasChanged(),
   getMinMeshWidth(),
   getMaxMeshWidth(),
   getNumberOfInnerVertices(),
   getNumberOfBoundaryVertices(),
   getNumberOfOuterVertices(),
   getNumberOfInnerCells(),
   getNumberOfOuterCells(),
   getNumberOfInnerLeafVertices(),
   getNumberOfBoundaryLeafVertices(),
   getNumberOfOuterLeafVertices(),
   getNumberOfInnerLeafCells(),
   getNumberOfOuterLeafCells(),
   getMaxLevel(),
   getHasRefined(),
   getHasTriggeredRefinementForNextIteration(),
   getHasErased(),
   getHasTriggeredEraseForNextIteration(),
   getHasChangedVertexOrCellState(),
   getHasModifiedGridInPreviousIteration(),
   getIsTraversalInverted(),
   getReduceStateAndCell(),
   getCouldNotEraseDueToDecompositionFlag(),
   getSubWorkerIsInvolvedInJoinOrFork()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );

MPI_Datatype exahype::records::State::Datatype = 0;
MPI_Datatype exahype::records::State::FullDatatype = 0;


void exahype::records::State::initDatatype() {
   {
      State dummyState[2];
      
      const int Attributes = 29;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //gridConstructionState
         MPI_CHAR,		 //firstGridSetupIteration
         MPI_INT,		 //mergeMode
         MPI_INT,		 //sendMode
         MPI_CHAR,		 //limiterDomainHasChanged
         MPI_DOUBLE,		 //minMeshWidth
         MPI_DOUBLE,		 //maxMeshWidth
         MPI_DOUBLE,		 //numberOfInnerVertices
         MPI_DOUBLE,		 //numberOfBoundaryVertices
         MPI_DOUBLE,		 //numberOfOuterVertices
         MPI_DOUBLE,		 //numberOfInnerCells
         MPI_DOUBLE,		 //numberOfOuterCells
         MPI_DOUBLE,		 //numberOfInnerLeafVertices
         MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
         MPI_DOUBLE,		 //numberOfOuterLeafVertices
         MPI_DOUBLE,		 //numberOfInnerLeafCells
         MPI_DOUBLE,		 //numberOfOuterLeafCells
         MPI_INT,		 //maxLevel
         MPI_CHAR,		 //hasRefined
         MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
         MPI_CHAR,		 //hasErased
         MPI_CHAR,		 //hasTriggeredEraseForNextIteration
         MPI_CHAR,		 //hasChangedVertexOrCellState
         MPI_CHAR,		 //hasModifiedGridInPreviousIteration
         MPI_CHAR,		 //isTraversalInverted
         MPI_CHAR,		 //reduceStateAndCell
         MPI_CHAR,		 //couldNotEraseDueToDecompositionFlag
         MPI_CHAR,		 //subWorkerIsInvolvedInJoinOrFork
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //gridConstructionState
         1,		 //firstGridSetupIteration
         1,		 //mergeMode
         1,		 //sendMode
         1,		 //limiterDomainHasChanged
         DIMENSIONS,		 //minMeshWidth
         DIMENSIONS,		 //maxMeshWidth
         1,		 //numberOfInnerVertices
         1,		 //numberOfBoundaryVertices
         1,		 //numberOfOuterVertices
         1,		 //numberOfInnerCells
         1,		 //numberOfOuterCells
         1,		 //numberOfInnerLeafVertices
         1,		 //numberOfBoundaryLeafVertices
         1,		 //numberOfOuterLeafVertices
         1,		 //numberOfInnerLeafCells
         1,		 //numberOfOuterLeafCells
         1,		 //maxLevel
         1,		 //hasRefined
         1,		 //hasTriggeredRefinementForNextIteration
         1,		 //hasErased
         1,		 //hasTriggeredEraseForNextIteration
         1,		 //hasChangedVertexOrCellState
         1,		 //hasModifiedGridInPreviousIteration
         1,		 //isTraversalInverted
         1,		 //reduceStateAndCell
         1,		 //couldNotEraseDueToDecompositionFlag
         1,		 //subWorkerIsInvolvedInJoinOrFork
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[4] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[5] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[6] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[7] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[8] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[9] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[10] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[11] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[12] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[13] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[14] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[15] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[16] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[17] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[18] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[19] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[20] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[21] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[22] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[23] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[24] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[25] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[26] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[27] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._gridConstructionState))), 		&disp[28] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype );
      MPI_Type_commit( &State::Datatype );
      
   }
   {
      State dummyState[2];
      
      const int Attributes = 32;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //gridConstructionState
         MPI_CHAR,		 //firstGridSetupIteration
         MPI_INT,		 //mergeMode
         MPI_INT,		 //sendMode
         MPI_CHAR,		 //reinitTimeStepData
         MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
         MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
         MPI_CHAR,		 //limiterDomainHasChanged
         MPI_DOUBLE,		 //minMeshWidth
         MPI_DOUBLE,		 //maxMeshWidth
         MPI_DOUBLE,		 //numberOfInnerVertices
         MPI_DOUBLE,		 //numberOfBoundaryVertices
         MPI_DOUBLE,		 //numberOfOuterVertices
         MPI_DOUBLE,		 //numberOfInnerCells
         MPI_DOUBLE,		 //numberOfOuterCells
         MPI_DOUBLE,		 //numberOfInnerLeafVertices
         MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
         MPI_DOUBLE,		 //numberOfOuterLeafVertices
         MPI_DOUBLE,		 //numberOfInnerLeafCells
         MPI_DOUBLE,		 //numberOfOuterLeafCells
         MPI_INT,		 //maxLevel
         MPI_CHAR,		 //hasRefined
         MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
         MPI_CHAR,		 //hasErased
         MPI_CHAR,		 //hasTriggeredEraseForNextIteration
         MPI_CHAR,		 //hasChangedVertexOrCellState
         MPI_CHAR,		 //hasModifiedGridInPreviousIteration
         MPI_CHAR,		 //isTraversalInverted
         MPI_CHAR,		 //reduceStateAndCell
         MPI_CHAR,		 //couldNotEraseDueToDecompositionFlag
         MPI_CHAR,		 //subWorkerIsInvolvedInJoinOrFork
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //gridConstructionState
         1,		 //firstGridSetupIteration
         1,		 //mergeMode
         1,		 //sendMode
         1,		 //reinitTimeStepData
         1,		 //stabilityConditionOfOneSolverWasViolated
         1,		 //timeStepSizeWeightForPredictionRerun
         1,		 //limiterDomainHasChanged
         DIMENSIONS,		 //minMeshWidth
         DIMENSIONS,		 //maxMeshWidth
         1,		 //numberOfInnerVertices
         1,		 //numberOfBoundaryVertices
         1,		 //numberOfOuterVertices
         1,		 //numberOfInnerCells
         1,		 //numberOfOuterCells
         1,		 //numberOfInnerLeafVertices
         1,		 //numberOfBoundaryLeafVertices
         1,		 //numberOfOuterLeafVertices
         1,		 //numberOfInnerLeafCells
         1,		 //numberOfOuterLeafCells
         1,		 //maxLevel
         1,		 //hasRefined
         1,		 //hasTriggeredRefinementForNextIteration
         1,		 //hasErased
         1,		 //hasTriggeredEraseForNextIteration
         1,		 //hasChangedVertexOrCellState
         1,		 //hasModifiedGridInPreviousIteration
         1,		 //isTraversalInverted
         1,		 //reduceStateAndCell
         1,		 //couldNotEraseDueToDecompositionFlag
         1,		 //subWorkerIsInvolvedInJoinOrFork
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reinitTimeStepData))), 		&disp[4] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[5] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[6] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[7] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[8] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[9] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[10] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[11] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[12] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[13] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[14] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[15] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[16] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[17] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[18] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[19] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[20] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[21] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[22] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[23] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[24] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[25] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[26] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[27] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[28] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[29] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[30] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._gridConstructionState))), 		&disp[31] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype );
      MPI_Type_commit( &State::FullDatatype );
      
   }
   
}


void exahype::records::State::shutdownDatatype() {
   MPI_Type_free( &State::Datatype );
   MPI_Type_free( &State::FullDatatype );
   
}

void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   _senderDestinationRank = destination;
   
   if (communicateSleep<0) {
   
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message exahype::records::State "
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
      msg << "was not able to send message exahype::records::State "
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
         msg << "testing for finished send task for exahype::records::State "
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
         "exahype::records::State",
         "send(int)", destination,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::State",
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



void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

   MPI_Status  status;
   const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
   _senderDestinationRank = status.MPI_SOURCE;
   if ( result != MPI_SUCCESS ) {
      std::ostringstream msg;
      msg << "failed to start to receive exahype::records::State from node "
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
      msg << "failed to start to receive exahype::records::State from node "
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
         msg << "testing for finished receive task for exahype::records::State failed: "
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
         "exahype::records::State",
         "receive(int)", source,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::State",
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



bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::State::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::PersistentRecords::PersistentRecords(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_gridConstructionState(gridConstructionState),
_firstGridSetupIteration(firstGridSetupIteration),
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_minMeshWidth(minMeshWidth),
_maxMeshWidth(maxMeshWidth),
_numberOfInnerVertices(numberOfInnerVertices),
_numberOfBoundaryVertices(numberOfBoundaryVertices),
_numberOfOuterVertices(numberOfOuterVertices),
_numberOfInnerCells(numberOfInnerCells),
_numberOfOuterCells(numberOfOuterCells),
_numberOfInnerLeafVertices(numberOfInnerLeafVertices),
_numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
_numberOfOuterLeafVertices(numberOfOuterLeafVertices),
_numberOfInnerLeafCells(numberOfInnerLeafCells),
_numberOfOuterLeafCells(numberOfOuterLeafCells),
_maxLevel(maxLevel),
_isTraversalInverted(isTraversalInverted) {
setHasRefined(hasRefined);
setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
setHasErased(hasErased);
setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
setReduceStateAndCell(reduceStateAndCell);
setCouldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag);
setSubWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork);
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}

exahype::records::StatePacked::StatePacked() {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._gridConstructionState, persistentRecords._firstGridSetupIteration, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted, persistentRecords.getReduceStateAndCell(), persistentRecords.getCouldNotEraseDueToDecompositionFlag(), persistentRecords.getSubWorkerIsInvolvedInJoinOrFork()) {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_persistentRecords(gridConstructionState, firstGridSetupIteration, mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::~StatePacked() { }

std::string exahype::records::StatePacked::toString(const GridConstructionState& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getGridConstructionStateMapping() {
return exahype::records::State::getGridConstructionStateMapping();
}

std::string exahype::records::StatePacked::toString(const MergeMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getMergeModeMapping() {
return exahype::records::State::getMergeModeMapping();
}

std::string exahype::records::StatePacked::toString(const SendMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getSendModeMapping() {
return exahype::records::State::getSendModeMapping();
}



std::string exahype::records::StatePacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::StatePacked::toString (std::ostream& out) const {
out << "("; 
out << "gridConstructionState:" << toString(getGridConstructionState());
out << ",";
out << "firstGridSetupIteration:" << getFirstGridSetupIteration();
out << ",";
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
out << ",";
out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
out << ",";
out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
out << ",";
out << "numberOfInnerCells:" << getNumberOfInnerCells();
out << ",";
out << "numberOfOuterCells:" << getNumberOfOuterCells();
out << ",";
out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
out << ",";
out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
out << ",";
out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
out << ",";
out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
out << ",";
out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
out << ",";
out << "maxLevel:" << getMaxLevel();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out << ",";
out << "reduceStateAndCell:" << getReduceStateAndCell();
out << ",";
out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
out << ",";
out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
out <<  ")";
}


exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::State exahype::records::StatePacked::convert() const{
return State(
getGridConstructionState(),
getFirstGridSetupIteration(),
getMergeMode(),
getSendMode(),
getReinitTimeStepData(),
getStabilityConditionOfOneSolverWasViolated(),
getTimeStepSizeWeightForPredictionRerun(),
getLimiterDomainHasChanged(),
getMinMeshWidth(),
getMaxMeshWidth(),
getNumberOfInnerVertices(),
getNumberOfBoundaryVertices(),
getNumberOfOuterVertices(),
getNumberOfInnerCells(),
getNumberOfOuterCells(),
getNumberOfInnerLeafVertices(),
getNumberOfBoundaryLeafVertices(),
getNumberOfOuterLeafVertices(),
getNumberOfInnerLeafCells(),
getNumberOfOuterLeafCells(),
getMaxLevel(),
getHasRefined(),
getHasTriggeredRefinementForNextIteration(),
getHasErased(),
getHasTriggeredEraseForNextIteration(),
getHasChangedVertexOrCellState(),
getHasModifiedGridInPreviousIteration(),
getIsTraversalInverted(),
getReduceStateAndCell(),
getCouldNotEraseDueToDecompositionFlag(),
getSubWorkerIsInvolvedInJoinOrFork()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );

MPI_Datatype exahype::records::StatePacked::Datatype = 0;
MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;


void exahype::records::StatePacked::initDatatype() {
{
   StatePacked dummyStatePacked[2];
   
   const int Attributes = 21;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //gridConstructionState
      MPI_CHAR,		 //firstGridSetupIteration
      MPI_INT,		 //mergeMode
      MPI_INT,		 //sendMode
      MPI_CHAR,		 //limiterDomainHasChanged
      MPI_DOUBLE,		 //minMeshWidth
      MPI_DOUBLE,		 //maxMeshWidth
      MPI_DOUBLE,		 //numberOfInnerVertices
      MPI_DOUBLE,		 //numberOfBoundaryVertices
      MPI_DOUBLE,		 //numberOfOuterVertices
      MPI_DOUBLE,		 //numberOfInnerCells
      MPI_DOUBLE,		 //numberOfOuterCells
      MPI_DOUBLE,		 //numberOfInnerLeafVertices
      MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
      MPI_DOUBLE,		 //numberOfOuterLeafVertices
      MPI_DOUBLE,		 //numberOfInnerLeafCells
      MPI_DOUBLE,		 //numberOfOuterLeafCells
      MPI_INT,		 //maxLevel
      MPI_CHAR,		 //isTraversalInverted
      MPI_SHORT,		 //_packedRecords0
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //gridConstructionState
      1,		 //firstGridSetupIteration
      1,		 //mergeMode
      1,		 //sendMode
      1,		 //limiterDomainHasChanged
      DIMENSIONS,		 //minMeshWidth
      DIMENSIONS,		 //maxMeshWidth
      1,		 //numberOfInnerVertices
      1,		 //numberOfBoundaryVertices
      1,		 //numberOfOuterVertices
      1,		 //numberOfInnerCells
      1,		 //numberOfOuterCells
      1,		 //numberOfInnerLeafVertices
      1,		 //numberOfBoundaryLeafVertices
      1,		 //numberOfOuterLeafVertices
      1,		 //numberOfInnerLeafCells
      1,		 //numberOfOuterLeafCells
      1,		 //maxLevel
      1,		 //isTraversalInverted
      1,		 //_packedRecords0
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[3] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[4] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[5] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[6] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[7] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[8] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[9] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[10] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[11] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[12] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[13] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[14] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[15] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[16] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[17] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[18] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[19] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._gridConstructionState))), 		&disp[20] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype );
   MPI_Type_commit( &StatePacked::Datatype );
   
}
{
   StatePacked dummyStatePacked[2];
   
   const int Attributes = 24;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //gridConstructionState
      MPI_CHAR,		 //firstGridSetupIteration
      MPI_INT,		 //mergeMode
      MPI_INT,		 //sendMode
      MPI_CHAR,		 //reinitTimeStepData
      MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
      MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
      MPI_CHAR,		 //limiterDomainHasChanged
      MPI_DOUBLE,		 //minMeshWidth
      MPI_DOUBLE,		 //maxMeshWidth
      MPI_DOUBLE,		 //numberOfInnerVertices
      MPI_DOUBLE,		 //numberOfBoundaryVertices
      MPI_DOUBLE,		 //numberOfOuterVertices
      MPI_DOUBLE,		 //numberOfInnerCells
      MPI_DOUBLE,		 //numberOfOuterCells
      MPI_DOUBLE,		 //numberOfInnerLeafVertices
      MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
      MPI_DOUBLE,		 //numberOfOuterLeafVertices
      MPI_DOUBLE,		 //numberOfInnerLeafCells
      MPI_DOUBLE,		 //numberOfOuterLeafCells
      MPI_INT,		 //maxLevel
      MPI_CHAR,		 //isTraversalInverted
      MPI_SHORT,		 //_packedRecords0
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //gridConstructionState
      1,		 //firstGridSetupIteration
      1,		 //mergeMode
      1,		 //sendMode
      1,		 //reinitTimeStepData
      1,		 //stabilityConditionOfOneSolverWasViolated
      1,		 //timeStepSizeWeightForPredictionRerun
      1,		 //limiterDomainHasChanged
      DIMENSIONS,		 //minMeshWidth
      DIMENSIONS,		 //maxMeshWidth
      1,		 //numberOfInnerVertices
      1,		 //numberOfBoundaryVertices
      1,		 //numberOfOuterVertices
      1,		 //numberOfInnerCells
      1,		 //numberOfOuterCells
      1,		 //numberOfInnerLeafVertices
      1,		 //numberOfBoundaryLeafVertices
      1,		 //numberOfOuterLeafVertices
      1,		 //numberOfInnerLeafCells
      1,		 //numberOfOuterLeafCells
      1,		 //maxLevel
      1,		 //isTraversalInverted
      1,		 //_packedRecords0
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[3] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._reinitTimeStepData))), 		&disp[4] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[5] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[6] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[7] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[8] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[9] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[10] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[11] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[12] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[13] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[14] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[15] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[16] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[17] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[18] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[19] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[20] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[21] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[22] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._gridConstructionState))), 		&disp[23] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype );
   MPI_Type_commit( &StatePacked::FullDatatype );
   
}

}


void exahype::records::StatePacked::shutdownDatatype() {
MPI_Type_free( &StatePacked::Datatype );
MPI_Type_free( &StatePacked::FullDatatype );

}

void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

   const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
   if  (result!=MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message exahype::records::StatePacked "
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
   msg << "was not able to send message exahype::records::StatePacked "
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
      msg << "testing for finished send task for exahype::records::StatePacked "
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
      "exahype::records::StatePacked",
      "send(int)", destination,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "exahype::records::StatePacked",
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



void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
   std::ostringstream msg;
   msg << "failed to start to receive exahype::records::StatePacked from node "
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
   msg << "failed to start to receive exahype::records::StatePacked from node "
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
      msg << "testing for finished receive task for exahype::records::StatePacked failed: "
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
      "exahype::records::StatePacked",
      "receive(int)", source,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "exahype::records::StatePacked",
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



bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::StatePacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(TrackGridStatistics) && !defined(Parallel)
exahype::records::State::PersistentRecords::PersistentRecords() {

}


exahype::records::State::PersistentRecords::PersistentRecords(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_minMeshWidth(minMeshWidth),
_maxMeshWidth(maxMeshWidth),
_numberOfInnerVertices(numberOfInnerVertices),
_numberOfBoundaryVertices(numberOfBoundaryVertices),
_numberOfOuterVertices(numberOfOuterVertices),
_numberOfInnerCells(numberOfInnerCells),
_numberOfOuterCells(numberOfOuterCells),
_numberOfInnerLeafVertices(numberOfInnerLeafVertices),
_numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
_numberOfOuterLeafVertices(numberOfOuterLeafVertices),
_numberOfInnerLeafCells(numberOfInnerLeafCells),
_numberOfOuterLeafCells(numberOfOuterLeafCells),
_maxLevel(maxLevel),
_hasRefined(hasRefined),
_hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
_hasErased(hasErased),
_hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
_hasChangedVertexOrCellState(hasChangedVertexOrCellState),
_hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
_isTraversalInverted(isTraversalInverted) {

}

exahype::records::State::State() {

}


exahype::records::State::State(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted) {

}


exahype::records::State::State(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_persistentRecords(mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {

}


exahype::records::State::~State() { }

std::string exahype::records::State::toString(const MergeMode& param) {
switch (param) {
case MergeNothing: return "MergeNothing";
case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
case MergeFaceData: return "MergeFaceData";
case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
}
return "undefined";
}

std::string exahype::records::State::getMergeModeMapping() {
return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,BroadcastAndMergeTimeStepDataAndMergeFaceData=3)";
}
std::string exahype::records::State::toString(const SendMode& param) {
switch (param) {
case SendNothing: return "SendNothing";
case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
case SendFaceData: return "SendFaceData";
case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
}
return "undefined";
}

std::string exahype::records::State::getSendModeMapping() {
return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
}


std::string exahype::records::State::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::State::toString (std::ostream& out) const {
out << "("; 
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
out << ",";
out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
out << ",";
out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
out << ",";
out << "numberOfInnerCells:" << getNumberOfInnerCells();
out << ",";
out << "numberOfOuterCells:" << getNumberOfOuterCells();
out << ",";
out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
out << ",";
out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
out << ",";
out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
out << ",";
out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
out << ",";
out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
out << ",";
out << "maxLevel:" << getMaxLevel();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out <<  ")";
}


exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::StatePacked exahype::records::State::convert() const{
return StatePacked(
getMergeMode(),
getSendMode(),
getReinitTimeStepData(),
getStabilityConditionOfOneSolverWasViolated(),
getTimeStepSizeWeightForPredictionRerun(),
getLimiterDomainHasChanged(),
getMinMeshWidth(),
getMaxMeshWidth(),
getNumberOfInnerVertices(),
getNumberOfBoundaryVertices(),
getNumberOfOuterVertices(),
getNumberOfInnerCells(),
getNumberOfOuterCells(),
getNumberOfInnerLeafVertices(),
getNumberOfBoundaryLeafVertices(),
getNumberOfOuterLeafVertices(),
getNumberOfInnerLeafCells(),
getNumberOfOuterLeafCells(),
getMaxLevel(),
getHasRefined(),
getHasTriggeredRefinementForNextIteration(),
getHasErased(),
getHasTriggeredEraseForNextIteration(),
getHasChangedVertexOrCellState(),
getHasModifiedGridInPreviousIteration(),
getIsTraversalInverted()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );

MPI_Datatype exahype::records::State::Datatype = 0;
MPI_Datatype exahype::records::State::FullDatatype = 0;


void exahype::records::State::initDatatype() {
{
State dummyState[2];

const int Attributes = 24;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //limiterDomainHasChanged
MPI_DOUBLE,		 //minMeshWidth
MPI_DOUBLE,		 //maxMeshWidth
MPI_DOUBLE,		 //numberOfInnerVertices
MPI_DOUBLE,		 //numberOfBoundaryVertices
MPI_DOUBLE,		 //numberOfOuterVertices
MPI_DOUBLE,		 //numberOfInnerCells
MPI_DOUBLE,		 //numberOfOuterCells
MPI_DOUBLE,		 //numberOfInnerLeafVertices
MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
MPI_DOUBLE,		 //numberOfOuterLeafVertices
MPI_DOUBLE,		 //numberOfInnerLeafCells
MPI_DOUBLE,		 //numberOfOuterLeafCells
MPI_INT,		 //maxLevel
MPI_CHAR,		 //hasRefined
MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
MPI_CHAR,		 //hasErased
MPI_CHAR,		 //hasTriggeredEraseForNextIteration
MPI_CHAR,		 //hasChangedVertexOrCellState
MPI_CHAR,		 //hasModifiedGridInPreviousIteration
MPI_CHAR,		 //isTraversalInverted
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //mergeMode
1,		 //sendMode
1,		 //limiterDomainHasChanged
DIMENSIONS,		 //minMeshWidth
DIMENSIONS,		 //maxMeshWidth
1,		 //numberOfInnerVertices
1,		 //numberOfBoundaryVertices
1,		 //numberOfOuterVertices
1,		 //numberOfInnerCells
1,		 //numberOfOuterCells
1,		 //numberOfInnerLeafVertices
1,		 //numberOfBoundaryLeafVertices
1,		 //numberOfOuterLeafVertices
1,		 //numberOfInnerLeafCells
1,		 //numberOfOuterLeafCells
1,		 //maxLevel
1,		 //hasRefined
1,		 //hasTriggeredRefinementForNextIteration
1,		 //hasErased
1,		 //hasTriggeredEraseForNextIteration
1,		 //hasChangedVertexOrCellState
1,		 //hasModifiedGridInPreviousIteration
1,		 //isTraversalInverted
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[18] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[19] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[20] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[21] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[22] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._mergeMode))), 		&disp[23] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype );
MPI_Type_commit( &State::Datatype );

}
{
State dummyState[2];

const int Attributes = 27;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //reinitTimeStepData
MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
MPI_CHAR,		 //limiterDomainHasChanged
MPI_DOUBLE,		 //minMeshWidth
MPI_DOUBLE,		 //maxMeshWidth
MPI_DOUBLE,		 //numberOfInnerVertices
MPI_DOUBLE,		 //numberOfBoundaryVertices
MPI_DOUBLE,		 //numberOfOuterVertices
MPI_DOUBLE,		 //numberOfInnerCells
MPI_DOUBLE,		 //numberOfOuterCells
MPI_DOUBLE,		 //numberOfInnerLeafVertices
MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
MPI_DOUBLE,		 //numberOfOuterLeafVertices
MPI_DOUBLE,		 //numberOfInnerLeafCells
MPI_DOUBLE,		 //numberOfOuterLeafCells
MPI_INT,		 //maxLevel
MPI_CHAR,		 //hasRefined
MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
MPI_CHAR,		 //hasErased
MPI_CHAR,		 //hasTriggeredEraseForNextIteration
MPI_CHAR,		 //hasChangedVertexOrCellState
MPI_CHAR,		 //hasModifiedGridInPreviousIteration
MPI_CHAR,		 //isTraversalInverted
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //mergeMode
1,		 //sendMode
1,		 //reinitTimeStepData
1,		 //stabilityConditionOfOneSolverWasViolated
1,		 //timeStepSizeWeightForPredictionRerun
1,		 //limiterDomainHasChanged
DIMENSIONS,		 //minMeshWidth
DIMENSIONS,		 //maxMeshWidth
1,		 //numberOfInnerVertices
1,		 //numberOfBoundaryVertices
1,		 //numberOfOuterVertices
1,		 //numberOfInnerCells
1,		 //numberOfOuterCells
1,		 //numberOfInnerLeafVertices
1,		 //numberOfBoundaryLeafVertices
1,		 //numberOfOuterLeafVertices
1,		 //numberOfInnerLeafCells
1,		 //numberOfOuterLeafCells
1,		 //maxLevel
1,		 //hasRefined
1,		 //hasTriggeredRefinementForNextIteration
1,		 //hasErased
1,		 //hasTriggeredEraseForNextIteration
1,		 //hasChangedVertexOrCellState
1,		 //hasModifiedGridInPreviousIteration
1,		 //isTraversalInverted
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reinitTimeStepData))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[18] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[19] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[20] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[21] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[22] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[23] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[24] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[25] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._mergeMode))), 		&disp[26] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype );
MPI_Type_commit( &State::FullDatatype );

}

}


void exahype::records::State::shutdownDatatype() {
MPI_Type_free( &State::Datatype );
MPI_Type_free( &State::FullDatatype );

}

void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::State "
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
msg << "was not able to send message exahype::records::State "
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
msg << "testing for finished send task for exahype::records::State "
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
"exahype::records::State",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::State",
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



void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::State from node "
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
msg << "failed to start to receive exahype::records::State from node "
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
msg << "testing for finished receive task for exahype::records::State failed: "
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
"exahype::records::State",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::State",
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



bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::State::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
if ((6 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((6 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::PersistentRecords::PersistentRecords(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_minMeshWidth(minMeshWidth),
_maxMeshWidth(maxMeshWidth),
_numberOfInnerVertices(numberOfInnerVertices),
_numberOfBoundaryVertices(numberOfBoundaryVertices),
_numberOfOuterVertices(numberOfOuterVertices),
_numberOfInnerCells(numberOfInnerCells),
_numberOfOuterCells(numberOfOuterCells),
_numberOfInnerLeafVertices(numberOfInnerLeafVertices),
_numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
_numberOfOuterLeafVertices(numberOfOuterLeafVertices),
_numberOfInnerLeafCells(numberOfInnerLeafCells),
_numberOfOuterLeafCells(numberOfOuterLeafCells),
_maxLevel(maxLevel),
_isTraversalInverted(isTraversalInverted) {
setHasRefined(hasRefined);
setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
setHasErased(hasErased);
setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
if ((6 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((6 < (8 * sizeof(short int))));

}

exahype::records::StatePacked::StatePacked() {
if ((6 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((6 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted) {
if ((6 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((6 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
_persistentRecords(mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
if ((6 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((6 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::~StatePacked() { }

std::string exahype::records::StatePacked::toString(const MergeMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getMergeModeMapping() {
return exahype::records::State::getMergeModeMapping();
}

std::string exahype::records::StatePacked::toString(const SendMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getSendModeMapping() {
return exahype::records::State::getSendModeMapping();
}



std::string exahype::records::StatePacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::StatePacked::toString (std::ostream& out) const {
out << "("; 
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
out << ",";
out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
out << ",";
out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
out << ",";
out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
out << ",";
out << "numberOfInnerCells:" << getNumberOfInnerCells();
out << ",";
out << "numberOfOuterCells:" << getNumberOfOuterCells();
out << ",";
out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
out << ",";
out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
out << ",";
out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
out << ",";
out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
out << ",";
out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
out << ",";
out << "maxLevel:" << getMaxLevel();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out <<  ")";
}


exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::State exahype::records::StatePacked::convert() const{
return State(
getMergeMode(),
getSendMode(),
getReinitTimeStepData(),
getStabilityConditionOfOneSolverWasViolated(),
getTimeStepSizeWeightForPredictionRerun(),
getLimiterDomainHasChanged(),
getMinMeshWidth(),
getMaxMeshWidth(),
getNumberOfInnerVertices(),
getNumberOfBoundaryVertices(),
getNumberOfOuterVertices(),
getNumberOfInnerCells(),
getNumberOfOuterCells(),
getNumberOfInnerLeafVertices(),
getNumberOfBoundaryLeafVertices(),
getNumberOfOuterLeafVertices(),
getNumberOfInnerLeafCells(),
getNumberOfOuterLeafCells(),
getMaxLevel(),
getHasRefined(),
getHasTriggeredRefinementForNextIteration(),
getHasErased(),
getHasTriggeredEraseForNextIteration(),
getHasChangedVertexOrCellState(),
getHasModifiedGridInPreviousIteration(),
getIsTraversalInverted()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );

MPI_Datatype exahype::records::StatePacked::Datatype = 0;
MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;


void exahype::records::StatePacked::initDatatype() {
{
StatePacked dummyStatePacked[2];

const int Attributes = 19;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //limiterDomainHasChanged
MPI_DOUBLE,		 //minMeshWidth
MPI_DOUBLE,		 //maxMeshWidth
MPI_DOUBLE,		 //numberOfInnerVertices
MPI_DOUBLE,		 //numberOfBoundaryVertices
MPI_DOUBLE,		 //numberOfOuterVertices
MPI_DOUBLE,		 //numberOfInnerCells
MPI_DOUBLE,		 //numberOfOuterCells
MPI_DOUBLE,		 //numberOfInnerLeafVertices
MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
MPI_DOUBLE,		 //numberOfOuterLeafVertices
MPI_DOUBLE,		 //numberOfInnerLeafCells
MPI_DOUBLE,		 //numberOfOuterLeafCells
MPI_INT,		 //maxLevel
MPI_CHAR,		 //isTraversalInverted
MPI_SHORT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //mergeMode
1,		 //sendMode
1,		 //limiterDomainHasChanged
DIMENSIONS,		 //minMeshWidth
DIMENSIONS,		 //maxMeshWidth
1,		 //numberOfInnerVertices
1,		 //numberOfBoundaryVertices
1,		 //numberOfOuterVertices
1,		 //numberOfInnerCells
1,		 //numberOfOuterCells
1,		 //numberOfInnerLeafVertices
1,		 //numberOfBoundaryLeafVertices
1,		 //numberOfOuterLeafVertices
1,		 //numberOfInnerLeafCells
1,		 //numberOfOuterLeafCells
1,		 //maxLevel
1,		 //isTraversalInverted
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._mergeMode))), 		&disp[18] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype );
MPI_Type_commit( &StatePacked::Datatype );

}
{
StatePacked dummyStatePacked[2];

const int Attributes = 22;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //reinitTimeStepData
MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
MPI_CHAR,		 //limiterDomainHasChanged
MPI_DOUBLE,		 //minMeshWidth
MPI_DOUBLE,		 //maxMeshWidth
MPI_DOUBLE,		 //numberOfInnerVertices
MPI_DOUBLE,		 //numberOfBoundaryVertices
MPI_DOUBLE,		 //numberOfOuterVertices
MPI_DOUBLE,		 //numberOfInnerCells
MPI_DOUBLE,		 //numberOfOuterCells
MPI_DOUBLE,		 //numberOfInnerLeafVertices
MPI_DOUBLE,		 //numberOfBoundaryLeafVertices
MPI_DOUBLE,		 //numberOfOuterLeafVertices
MPI_DOUBLE,		 //numberOfInnerLeafCells
MPI_DOUBLE,		 //numberOfOuterLeafCells
MPI_INT,		 //maxLevel
MPI_CHAR,		 //isTraversalInverted
MPI_SHORT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //mergeMode
1,		 //sendMode
1,		 //reinitTimeStepData
1,		 //stabilityConditionOfOneSolverWasViolated
1,		 //timeStepSizeWeightForPredictionRerun
1,		 //limiterDomainHasChanged
DIMENSIONS,		 //minMeshWidth
DIMENSIONS,		 //maxMeshWidth
1,		 //numberOfInnerVertices
1,		 //numberOfBoundaryVertices
1,		 //numberOfOuterVertices
1,		 //numberOfInnerCells
1,		 //numberOfOuterCells
1,		 //numberOfInnerLeafVertices
1,		 //numberOfBoundaryLeafVertices
1,		 //numberOfOuterLeafVertices
1,		 //numberOfInnerLeafCells
1,		 //numberOfOuterLeafCells
1,		 //maxLevel
1,		 //isTraversalInverted
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._reinitTimeStepData))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[18] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[19] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[20] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._mergeMode))), 		&disp[21] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype );
MPI_Type_commit( &StatePacked::FullDatatype );

}

}


void exahype::records::StatePacked::shutdownDatatype() {
MPI_Type_free( &StatePacked::Datatype );
MPI_Type_free( &StatePacked::FullDatatype );

}

void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::StatePacked "
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
msg << "was not able to send message exahype::records::StatePacked "
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
msg << "testing for finished send task for exahype::records::StatePacked "
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
"exahype::records::StatePacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::StatePacked",
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



void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::StatePacked from node "
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
msg << "failed to start to receive exahype::records::StatePacked from node "
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
msg << "testing for finished receive task for exahype::records::StatePacked failed: "
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
"exahype::records::StatePacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::StatePacked",
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



bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::StatePacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif !defined(TrackGridStatistics) && defined(Parallel)
exahype::records::State::PersistentRecords::PersistentRecords() {

}


exahype::records::State::PersistentRecords::PersistentRecords(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_gridConstructionState(gridConstructionState),
_firstGridSetupIteration(firstGridSetupIteration),
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_hasRefined(hasRefined),
_hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
_hasErased(hasErased),
_hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
_hasChangedVertexOrCellState(hasChangedVertexOrCellState),
_hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
_isTraversalInverted(isTraversalInverted),
_reduceStateAndCell(reduceStateAndCell),
_couldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag),
_subWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork) {

}

exahype::records::State::State() {

}


exahype::records::State::State(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._gridConstructionState, persistentRecords._firstGridSetupIteration, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted, persistentRecords._reduceStateAndCell, persistentRecords._couldNotEraseDueToDecompositionFlag, persistentRecords._subWorkerIsInvolvedInJoinOrFork) {

}


exahype::records::State::State(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_persistentRecords(gridConstructionState, firstGridSetupIteration, mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {

}


exahype::records::State::~State() { }

std::string exahype::records::State::toString(const GridConstructionState& param) {
switch (param) {
case Default: return "Default";
case Veto: return "Veto";
case Aggressive: return "Aggressive";
}
return "undefined";
}

std::string exahype::records::State::getGridConstructionStateMapping() {
return "GridConstructionState(Default=0,Veto=1,Aggressive=2)";
}
std::string exahype::records::State::toString(const MergeMode& param) {
switch (param) {
case MergeNothing: return "MergeNothing";
case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
case MergeFaceData: return "MergeFaceData";
case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
}
return "undefined";
}

std::string exahype::records::State::getMergeModeMapping() {
return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,BroadcastAndMergeTimeStepDataAndMergeFaceData=3)";
}
std::string exahype::records::State::toString(const SendMode& param) {
switch (param) {
case SendNothing: return "SendNothing";
case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
case SendFaceData: return "SendFaceData";
case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
}
return "undefined";
}

std::string exahype::records::State::getSendModeMapping() {
return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
}


std::string exahype::records::State::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::State::toString (std::ostream& out) const {
out << "("; 
out << "gridConstructionState:" << toString(getGridConstructionState());
out << ",";
out << "firstGridSetupIteration:" << getFirstGridSetupIteration();
out << ",";
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out << ",";
out << "reduceStateAndCell:" << getReduceStateAndCell();
out << ",";
out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
out << ",";
out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
out <<  ")";
}


exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::StatePacked exahype::records::State::convert() const{
return StatePacked(
getGridConstructionState(),
getFirstGridSetupIteration(),
getMergeMode(),
getSendMode(),
getReinitTimeStepData(),
getStabilityConditionOfOneSolverWasViolated(),
getTimeStepSizeWeightForPredictionRerun(),
getLimiterDomainHasChanged(),
getHasRefined(),
getHasTriggeredRefinementForNextIteration(),
getHasErased(),
getHasTriggeredEraseForNextIteration(),
getHasChangedVertexOrCellState(),
getHasModifiedGridInPreviousIteration(),
getIsTraversalInverted(),
getReduceStateAndCell(),
getCouldNotEraseDueToDecompositionFlag(),
getSubWorkerIsInvolvedInJoinOrFork()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );

MPI_Datatype exahype::records::State::Datatype = 0;
MPI_Datatype exahype::records::State::FullDatatype = 0;


void exahype::records::State::initDatatype() {
{
State dummyState[2];

const int Attributes = 16;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //gridConstructionState
MPI_CHAR,		 //firstGridSetupIteration
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //limiterDomainHasChanged
MPI_CHAR,		 //hasRefined
MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
MPI_CHAR,		 //hasErased
MPI_CHAR,		 //hasTriggeredEraseForNextIteration
MPI_CHAR,		 //hasChangedVertexOrCellState
MPI_CHAR,		 //hasModifiedGridInPreviousIteration
MPI_CHAR,		 //isTraversalInverted
MPI_CHAR,		 //reduceStateAndCell
MPI_CHAR,		 //couldNotEraseDueToDecompositionFlag
MPI_CHAR,		 //subWorkerIsInvolvedInJoinOrFork
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //gridConstructionState
1,		 //firstGridSetupIteration
1,		 //mergeMode
1,		 //sendMode
1,		 //limiterDomainHasChanged
1,		 //hasRefined
1,		 //hasTriggeredRefinementForNextIteration
1,		 //hasErased
1,		 //hasTriggeredEraseForNextIteration
1,		 //hasChangedVertexOrCellState
1,		 //hasModifiedGridInPreviousIteration
1,		 //isTraversalInverted
1,		 //reduceStateAndCell
1,		 //couldNotEraseDueToDecompositionFlag
1,		 //subWorkerIsInvolvedInJoinOrFork
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._gridConstructionState))), 		&disp[15] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype );
MPI_Type_commit( &State::Datatype );

}
{
State dummyState[2];

const int Attributes = 19;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //gridConstructionState
MPI_CHAR,		 //firstGridSetupIteration
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //reinitTimeStepData
MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
MPI_CHAR,		 //limiterDomainHasChanged
MPI_CHAR,		 //hasRefined
MPI_CHAR,		 //hasTriggeredRefinementForNextIteration
MPI_CHAR,		 //hasErased
MPI_CHAR,		 //hasTriggeredEraseForNextIteration
MPI_CHAR,		 //hasChangedVertexOrCellState
MPI_CHAR,		 //hasModifiedGridInPreviousIteration
MPI_CHAR,		 //isTraversalInverted
MPI_CHAR,		 //reduceStateAndCell
MPI_CHAR,		 //couldNotEraseDueToDecompositionFlag
MPI_CHAR,		 //subWorkerIsInvolvedInJoinOrFork
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //gridConstructionState
1,		 //firstGridSetupIteration
1,		 //mergeMode
1,		 //sendMode
1,		 //reinitTimeStepData
1,		 //stabilityConditionOfOneSolverWasViolated
1,		 //timeStepSizeWeightForPredictionRerun
1,		 //limiterDomainHasChanged
1,		 //hasRefined
1,		 //hasTriggeredRefinementForNextIteration
1,		 //hasErased
1,		 //hasTriggeredEraseForNextIteration
1,		 //hasChangedVertexOrCellState
1,		 //hasModifiedGridInPreviousIteration
1,		 //isTraversalInverted
1,		 //reduceStateAndCell
1,		 //couldNotEraseDueToDecompositionFlag
1,		 //subWorkerIsInvolvedInJoinOrFork
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reinitTimeStepData))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]._persistentRecords._gridConstructionState))), 		&disp[18] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype );
MPI_Type_commit( &State::FullDatatype );

}

}


void exahype::records::State::shutdownDatatype() {
MPI_Type_free( &State::Datatype );
MPI_Type_free( &State::FullDatatype );

}

void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::State "
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
msg << "was not able to send message exahype::records::State "
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
msg << "testing for finished send task for exahype::records::State "
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
"exahype::records::State",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::State",
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



void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::State from node "
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
msg << "failed to start to receive exahype::records::State from node "
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
msg << "testing for finished receive task for exahype::records::State failed: "
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
"exahype::records::State",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::State",
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



bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::State::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::PersistentRecords::PersistentRecords(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_gridConstructionState(gridConstructionState),
_firstGridSetupIteration(firstGridSetupIteration),
_mergeMode(mergeMode),
_sendMode(sendMode),
_reinitTimeStepData(reinitTimeStepData),
_stabilityConditionOfOneSolverWasViolated(stabilityConditionOfOneSolverWasViolated),
_timeStepSizeWeightForPredictionRerun(timeStepSizeWeightForPredictionRerun),
_limiterDomainHasChanged(limiterDomainHasChanged),
_isTraversalInverted(isTraversalInverted) {
setHasRefined(hasRefined);
setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
setHasErased(hasErased);
setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
setReduceStateAndCell(reduceStateAndCell);
setCouldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag);
setSubWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork);
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}

exahype::records::StatePacked::StatePacked() {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._gridConstructionState, persistentRecords._firstGridSetupIteration, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._reinitTimeStepData, persistentRecords._stabilityConditionOfOneSolverWasViolated, persistentRecords._timeStepSizeWeightForPredictionRerun, persistentRecords._limiterDomainHasChanged, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted, persistentRecords.getReduceStateAndCell(), persistentRecords.getCouldNotEraseDueToDecompositionFlag(), persistentRecords.getSubWorkerIsInvolvedInJoinOrFork()) {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::StatePacked(const GridConstructionState& gridConstructionState, const bool& firstGridSetupIteration, const MergeMode& mergeMode, const SendMode& sendMode, const bool& reinitTimeStepData, const bool& stabilityConditionOfOneSolverWasViolated, const double& timeStepSizeWeightForPredictionRerun, const bool& limiterDomainHasChanged, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
_persistentRecords(gridConstructionState, firstGridSetupIteration, mergeMode, sendMode, reinitTimeStepData, stabilityConditionOfOneSolverWasViolated, timeStepSizeWeightForPredictionRerun, limiterDomainHasChanged, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
if ((9 >= (8 * sizeof(short int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((9 < (8 * sizeof(short int))));

}


exahype::records::StatePacked::~StatePacked() { }

std::string exahype::records::StatePacked::toString(const GridConstructionState& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getGridConstructionStateMapping() {
return exahype::records::State::getGridConstructionStateMapping();
}

std::string exahype::records::StatePacked::toString(const MergeMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getMergeModeMapping() {
return exahype::records::State::getMergeModeMapping();
}

std::string exahype::records::StatePacked::toString(const SendMode& param) {
return exahype::records::State::toString(param);
}

std::string exahype::records::StatePacked::getSendModeMapping() {
return exahype::records::State::getSendModeMapping();
}



std::string exahype::records::StatePacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::StatePacked::toString (std::ostream& out) const {
out << "("; 
out << "gridConstructionState:" << toString(getGridConstructionState());
out << ",";
out << "firstGridSetupIteration:" << getFirstGridSetupIteration();
out << ",";
out << "mergeMode:" << toString(getMergeMode());
out << ",";
out << "sendMode:" << toString(getSendMode());
out << ",";
out << "reinitTimeStepData:" << getReinitTimeStepData();
out << ",";
out << "stabilityConditionOfOneSolverWasViolated:" << getStabilityConditionOfOneSolverWasViolated();
out << ",";
out << "timeStepSizeWeightForPredictionRerun:" << getTimeStepSizeWeightForPredictionRerun();
out << ",";
out << "limiterDomainHasChanged:" << getLimiterDomainHasChanged();
out << ",";
out << "hasRefined:" << getHasRefined();
out << ",";
out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
out << ",";
out << "hasErased:" << getHasErased();
out << ",";
out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
out << ",";
out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
out << ",";
out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
out << ",";
out << "isTraversalInverted:" << getIsTraversalInverted();
out << ",";
out << "reduceStateAndCell:" << getReduceStateAndCell();
out << ",";
out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
out << ",";
out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
out <<  ")";
}


exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::State exahype::records::StatePacked::convert() const{
return State(
getGridConstructionState(),
getFirstGridSetupIteration(),
getMergeMode(),
getSendMode(),
getReinitTimeStepData(),
getStabilityConditionOfOneSolverWasViolated(),
getTimeStepSizeWeightForPredictionRerun(),
getLimiterDomainHasChanged(),
getHasRefined(),
getHasTriggeredRefinementForNextIteration(),
getHasErased(),
getHasTriggeredEraseForNextIteration(),
getHasChangedVertexOrCellState(),
getHasModifiedGridInPreviousIteration(),
getIsTraversalInverted(),
getReduceStateAndCell(),
getCouldNotEraseDueToDecompositionFlag(),
getSubWorkerIsInvolvedInJoinOrFork()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );

MPI_Datatype exahype::records::StatePacked::Datatype = 0;
MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;


void exahype::records::StatePacked::initDatatype() {
{
StatePacked dummyStatePacked[2];

const int Attributes = 8;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //gridConstructionState
MPI_CHAR,		 //firstGridSetupIteration
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //limiterDomainHasChanged
MPI_CHAR,		 //isTraversalInverted
MPI_SHORT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //gridConstructionState
1,		 //firstGridSetupIteration
1,		 //mergeMode
1,		 //sendMode
1,		 //limiterDomainHasChanged
1,		 //isTraversalInverted
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._gridConstructionState))), 		&disp[7] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype );
MPI_Type_commit( &StatePacked::Datatype );

}
{
StatePacked dummyStatePacked[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //gridConstructionState
MPI_CHAR,		 //firstGridSetupIteration
MPI_INT,		 //mergeMode
MPI_INT,		 //sendMode
MPI_CHAR,		 //reinitTimeStepData
MPI_CHAR,		 //stabilityConditionOfOneSolverWasViolated
MPI_DOUBLE,		 //timeStepSizeWeightForPredictionRerun
MPI_CHAR,		 //limiterDomainHasChanged
MPI_CHAR,		 //isTraversalInverted
MPI_SHORT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //gridConstructionState
1,		 //firstGridSetupIteration
1,		 //mergeMode
1,		 //sendMode
1,		 //reinitTimeStepData
1,		 //stabilityConditionOfOneSolverWasViolated
1,		 //timeStepSizeWeightForPredictionRerun
1,		 //limiterDomainHasChanged
1,		 //isTraversalInverted
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._gridConstructionState))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._firstGridSetupIteration))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._reinitTimeStepData))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._stabilityConditionOfOneSolverWasViolated))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._timeStepSizeWeightForPredictionRerun))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._limiterDomainHasChanged))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]._persistentRecords._gridConstructionState))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype );
MPI_Type_commit( &StatePacked::FullDatatype );

}

}


void exahype::records::StatePacked::shutdownDatatype() {
MPI_Type_free( &StatePacked::Datatype );
MPI_Type_free( &StatePacked::FullDatatype );

}

void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::StatePacked "
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
msg << "was not able to send message exahype::records::StatePacked "
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
msg << "testing for finished send task for exahype::records::StatePacked "
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
"exahype::records::StatePacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::StatePacked",
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



void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::StatePacked from node "
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
msg << "failed to start to receive exahype::records::StatePacked from node "
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
msg << "testing for finished receive task for exahype::records::StatePacked failed: "
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
"exahype::records::StatePacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::StatePacked",
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



bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::StatePacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#endif


