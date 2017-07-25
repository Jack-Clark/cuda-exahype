#include "exahype/records/ADERDGCellDescription.h"

#if defined(Parallel)
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForNeighbourCommunication, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _riemannSolvePerformed(riemannSolvePerformed),
   _isInside(isInside),
   _adjacentToRemoteRank(adjacentToRemoteRank),
   _hasToHoldDataForNeighbourCommunication(hasToHoldDataForNeighbourCommunication),
   _hasToHoldDataForMasterWorkerCommunication(hasToHoldDataForMasterWorkerCommunication),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _parentIndex(parentIndex),
   _type(type),
   _refinementEvent(refinementEvent),
   _level(level),
   _offset(offset),
   _size(size),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _limiterStatus(limiterStatus),
   _mergedLimiterStatus(mergedLimiterStatus),
   _compressionState(compressionState),
   _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
   _bytesPerDoFInSolution(bytesPerDoFInSolution),
   _bytesPerDoFInUpdate(bytesPerDoFInUpdate),
   _bytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor),
   _bytesPerDoFInFluctuation(bytesPerDoFInFluctuation) {
      
   }
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords._riemannSolvePerformed, persistentRecords._isInside, persistentRecords._adjacentToRemoteRank, persistentRecords._hasToHoldDataForNeighbourCommunication, persistentRecords._hasToHoldDataForMasterWorkerCommunication, persistentRecords._faceDataExchangeCounter, persistentRecords._parentIndex, persistentRecords._type, persistentRecords._refinementEvent, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._limiterStatus, persistentRecords._mergedLimiterStatus, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInUpdate, persistentRecords._bytesPerDoFInExtrapolatedPredictor, persistentRecords._bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForNeighbourCommunication, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, riemannSolvePerformed, isInside, adjacentToRemoteRank, hasToHoldDataForNeighbourCommunication, hasToHoldDataForMasterWorkerCommunication, faceDataExchangeCounter, parentIndex, type, refinementEvent, level, offset, size, correctorTimeStepSize, correctorTimeStamp, previousCorrectorTimeStepSize, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, limiterStatus, mergedLimiterStatus, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::~ADERDGCellDescription() { }
   
   std::string exahype::records::ADERDGCellDescription::toString(const CompressionState& param) {
      switch (param) {
         case Uncompressed: return "Uncompressed";
         case CurrentlyProcessed: return "CurrentlyProcessed";
         case Compressed: return "Compressed";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getCompressionStateMapping() {
      return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const LimiterStatus& param) {
      switch (param) {
         case Ok: return "Ok";
         case NeighbourIsNeighbourOfTroubledCell: return "NeighbourIsNeighbourOfTroubledCell";
         case NeighbourIsTroubledCell: return "NeighbourIsTroubledCell";
         case Troubled: return "Troubled";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getLimiterStatusMapping() {
      return "LimiterStatus(Ok=0,NeighbourIsNeighbourOfTroubledCell=1,NeighbourIsTroubledCell=2,Troubled=3)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const RefinementEvent& param) {
      switch (param) {
         case None: return "None";
         case ErasingChildrenRequested: return "ErasingChildrenRequested";
         case ErasingChildren: return "ErasingChildren";
         case ChangeChildrenToDescendantsRequested: return "ChangeChildrenToDescendantsRequested";
         case ChangeChildrenToDescendants: return "ChangeChildrenToDescendants";
         case RefiningRequested: return "RefiningRequested";
         case Refining: return "Refining";
         case DeaugmentingChildrenRequestedTriggered: return "DeaugmentingChildrenRequestedTriggered";
         case DeaugmentingChildrenRequested: return "DeaugmentingChildrenRequested";
         case DeaugmentingChildren: return "DeaugmentingChildren";
         case AugmentingRequested: return "AugmentingRequested";
         case Augmenting: return "Augmenting";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getRefinementEventMapping() {
      return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToDescendantsRequested=3,ChangeChildrenToDescendants=4,RefiningRequested=5,Refining=6,DeaugmentingChildrenRequestedTriggered=7,DeaugmentingChildrenRequested=8,DeaugmentingChildren=9,AugmentingRequested=10,Augmenting=11)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const Type& param) {
      switch (param) {
         case Erased: return "Erased";
         case Ancestor: return "Ancestor";
         case EmptyAncestor: return "EmptyAncestor";
         case Cell: return "Cell";
         case Descendant: return "Descendant";
         case EmptyDescendant: return "EmptyDescendant";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getTypeMapping() {
      return "Type(Erased=0,Ancestor=1,EmptyAncestor=2,Cell=3,Descendant=4,EmptyDescendant=5)";
   }
   
   
   std::string exahype::records::ADERDGCellDescription::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescription::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "riemannSolvePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getRiemannSolvePerformed(i) << ",";
   }
   out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "isInside:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getIsInside(i) << ",";
   }
   out << getIsInside(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "adjacentToRemoteRank:" << getAdjacentToRemoteRank();
      out << ",";
      out << "hasToHoldDataForNeighbourCommunication:" << getHasToHoldDataForNeighbourCommunication();
      out << ",";
      out << "hasToHoldDataForMasterWorkerCommunication:" << getHasToHoldDataForMasterWorkerCommunication();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "limiterStatus:" << toString(getLimiterStatus());
      out << ",";
      out << "mergedLimiterStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getMergedLimiterStatus(i) << ",";
   }
   out << getMergedLimiterStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords exahype::records::ADERDGCellDescription::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescriptionPacked exahype::records::ADERDGCellDescription::convert() const{
      return ADERDGCellDescriptionPacked(
         getSolverNumber(),
         getRiemannSolvePerformed(),
         getIsInside(),
         getAdjacentToRemoteRank(),
         getHasToHoldDataForNeighbourCommunication(),
         getHasToHoldDataForMasterWorkerCommunication(),
         getFaceDataExchangeCounter(),
         getParentIndex(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getLimiterStatus(),
         getMergedLimiterStatus(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescription::_log( "exahype::records::ADERDGCellDescription" );
      
      MPI_Datatype exahype::records::ADERDGCellDescription::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescription::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescription::initDatatype() {
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            const int Attributes = 44;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //riemannSolvePerformed
               MPI_INT,		 //isInside
               MPI_CHAR,		 //adjacentToRemoteRank
               MPI_CHAR,		 //hasToHoldDataForNeighbourCommunication
               MPI_CHAR,		 //hasToHoldDataForMasterWorkerCommunication
               MPI_INT,		 //faceDataExchangeCounter
               MPI_INT,		 //parentIndex
               MPI_INT,		 //type
               MPI_INT,		 //refinementEvent
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //compressionState
               MPI_INT,		 //bytesPerDoFInPreviousSolution
               MPI_INT,		 //bytesPerDoFInSolution
               MPI_INT,		 //bytesPerDoFInUpdate
               MPI_INT,		 //bytesPerDoFInExtrapolatedPredictor
               MPI_INT,		 //bytesPerDoFInFluctuation
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               DIMENSIONS_TIMES_TWO,		 //riemannSolvePerformed
               DIMENSIONS_TIMES_TWO,		 //isInside
               1,		 //adjacentToRemoteRank
               1,		 //hasToHoldDataForNeighbourCommunication
               1,		 //hasToHoldDataForMasterWorkerCommunication
               DIMENSIONS_TIMES_TWO,		 //faceDataExchangeCounter
               1,		 //parentIndex
               1,		 //type
               1,		 //refinementEvent
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //compressionState
               1,		 //bytesPerDoFInPreviousSolution
               1,		 //bytesPerDoFInSolution
               1,		 //bytesPerDoFInUpdate
               1,		 //bytesPerDoFInExtrapolatedPredictor
               1,		 //bytesPerDoFInFluctuation
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._riemannSolvePerformed))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._isInside))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._adjacentToRemoteRank))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForNeighbourCommunication))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[32] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[33] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[34] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._limiterStatus))), 		&disp[35] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[36] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[37] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[38] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[39] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[40] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[41] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[42] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))), 		&disp[43] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::Datatype );
            MPI_Type_commit( &ADERDGCellDescription::Datatype );
            
         }
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            const int Attributes = 44;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //riemannSolvePerformed
               MPI_INT,		 //isInside
               MPI_CHAR,		 //adjacentToRemoteRank
               MPI_CHAR,		 //hasToHoldDataForNeighbourCommunication
               MPI_CHAR,		 //hasToHoldDataForMasterWorkerCommunication
               MPI_INT,		 //faceDataExchangeCounter
               MPI_INT,		 //parentIndex
               MPI_INT,		 //type
               MPI_INT,		 //refinementEvent
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //compressionState
               MPI_INT,		 //bytesPerDoFInPreviousSolution
               MPI_INT,		 //bytesPerDoFInSolution
               MPI_INT,		 //bytesPerDoFInUpdate
               MPI_INT,		 //bytesPerDoFInExtrapolatedPredictor
               MPI_INT,		 //bytesPerDoFInFluctuation
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               DIMENSIONS_TIMES_TWO,		 //riemannSolvePerformed
               DIMENSIONS_TIMES_TWO,		 //isInside
               1,		 //adjacentToRemoteRank
               1,		 //hasToHoldDataForNeighbourCommunication
               1,		 //hasToHoldDataForMasterWorkerCommunication
               DIMENSIONS_TIMES_TWO,		 //faceDataExchangeCounter
               1,		 //parentIndex
               1,		 //type
               1,		 //refinementEvent
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //compressionState
               1,		 //bytesPerDoFInPreviousSolution
               1,		 //bytesPerDoFInSolution
               1,		 //bytesPerDoFInUpdate
               1,		 //bytesPerDoFInExtrapolatedPredictor
               1,		 //bytesPerDoFInFluctuation
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._riemannSolvePerformed))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._isInside))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._adjacentToRemoteRank))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForNeighbourCommunication))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[32] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[33] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[34] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._limiterStatus))), 		&disp[35] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[36] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[37] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[38] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[39] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[40] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[41] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[42] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))), 		&disp[43] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescription::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescription::Datatype );
         MPI_Type_free( &ADERDGCellDescription::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::ADERDGCellDescription "
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
               msg << "was not able to send message exahype::records::ADERDGCellDescription "
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
                  msg << "testing for finished send task for exahype::records::ADERDGCellDescription "
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
                  "exahype::records::ADERDGCellDescription",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescription",
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
      
      
      
      void exahype::records::ADERDGCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::ADERDGCellDescription from node "
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
               msg << "failed to start to receive exahype::records::ADERDGCellDescription from node "
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
                  msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: "
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
                  "exahype::records::ADERDGCellDescription",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescription",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::ADERDGCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      
      
   #endif
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords() {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForNeighbourCommunication, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _hasToHoldDataForMasterWorkerCommunication(hasToHoldDataForMasterWorkerCommunication),
   _faceDataExchangeCounter(faceDataExchangeCounter),
   _parentIndex(parentIndex),
   _level(level),
   _offset(offset),
   _size(size),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _limiterStatus(limiterStatus),
   _mergedLimiterStatus(mergedLimiterStatus) {
      setRiemannSolvePerformed(riemannSolvePerformed);
      setIsInside(isInside);
      setAdjacentToRemoteRank(adjacentToRemoteRank);
      setHasToHoldDataForNeighbourCommunication(hasToHoldDataForNeighbourCommunication);
      setType(type);
      setRefinementEvent(refinementEvent);
      setCompressionState(compressionState);
      setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
      setBytesPerDoFInSolution(bytesPerDoFInSolution);
      setBytesPerDoFInUpdate(bytesPerDoFInUpdate);
      setBytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor);
      setBytesPerDoFInFluctuation(bytesPerDoFInFluctuation);
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 < (8 * sizeof(int))));
      
   }
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords.getRiemannSolvePerformed(), persistentRecords.getIsInside(), persistentRecords.getAdjacentToRemoteRank(), persistentRecords.getHasToHoldDataForNeighbourCommunication(), persistentRecords._hasToHoldDataForMasterWorkerCommunication, persistentRecords._faceDataExchangeCounter, persistentRecords._parentIndex, persistentRecords.getType(), persistentRecords.getRefinementEvent(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._limiterStatus, persistentRecords._mergedLimiterStatus, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInUpdate(), persistentRecords.getBytesPerDoFInExtrapolatedPredictor(), persistentRecords.getBytesPerDoFInFluctuation()) {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const bool& adjacentToRemoteRank, const bool& hasToHoldDataForNeighbourCommunication, const bool& hasToHoldDataForMasterWorkerCommunication, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,int>& faceDataExchangeCounter, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, riemannSolvePerformed, isInside, adjacentToRemoteRank, hasToHoldDataForNeighbourCommunication, hasToHoldDataForMasterWorkerCommunication, faceDataExchangeCounter, parentIndex, type, refinementEvent, level, offset, size, correctorTimeStepSize, correctorTimeStamp, previousCorrectorTimeStepSize, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, limiterStatus, mergedLimiterStatus, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+26 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() { }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Type& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getTypeMapping() {
      return exahype::records::ADERDGCellDescription::getTypeMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const RefinementEvent& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getRefinementEventMapping() {
      return exahype::records::ADERDGCellDescription::getRefinementEventMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const LimiterStatus& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getLimiterStatusMapping() {
      return exahype::records::ADERDGCellDescription::getLimiterStatusMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const CompressionState& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getCompressionStateMapping() {
      return exahype::records::ADERDGCellDescription::getCompressionStateMapping();
   }
   
   
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescriptionPacked::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "riemannSolvePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getRiemannSolvePerformed(i) << ",";
   }
   out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "isInside:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getIsInside(i) << ",";
   }
   out << getIsInside(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "adjacentToRemoteRank:" << getAdjacentToRemoteRank();
      out << ",";
      out << "hasToHoldDataForNeighbourCommunication:" << getHasToHoldDataForNeighbourCommunication();
      out << ",";
      out << "hasToHoldDataForMasterWorkerCommunication:" << getHasToHoldDataForMasterWorkerCommunication();
      out << ",";
      out << "faceDataExchangeCounter:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getFaceDataExchangeCounter(i) << ",";
   }
   out << getFaceDataExchangeCounter(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "limiterStatus:" << toString(getLimiterStatus());
      out << ",";
      out << "mergedLimiterStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getMergedLimiterStatus(i) << ",";
   }
   out << getMergedLimiterStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords exahype::records::ADERDGCellDescriptionPacked::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescription exahype::records::ADERDGCellDescriptionPacked::convert() const{
      return ADERDGCellDescription(
         getSolverNumber(),
         getRiemannSolvePerformed(),
         getIsInside(),
         getAdjacentToRemoteRank(),
         getHasToHoldDataForNeighbourCommunication(),
         getHasToHoldDataForMasterWorkerCommunication(),
         getFaceDataExchangeCounter(),
         getParentIndex(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getLimiterStatus(),
         getMergedLimiterStatus(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescriptionPacked::_log( "exahype::records::ADERDGCellDescriptionPacked" );
      
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescriptionPacked::initDatatype() {
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            const int Attributes = 33;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_CHAR,		 //hasToHoldDataForMasterWorkerCommunication
               MPI_INT,		 //faceDataExchangeCounter
               MPI_INT,		 //parentIndex
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //_packedRecords0
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               1,		 //hasToHoldDataForMasterWorkerCommunication
               DIMENSIONS_TIMES_TWO,		 //faceDataExchangeCounter
               1,		 //parentIndex
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //_packedRecords0
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._limiterStatus))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[32] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::Datatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
            
         }
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            const int Attributes = 33;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_CHAR,		 //hasToHoldDataForMasterWorkerCommunication
               MPI_INT,		 //faceDataExchangeCounter
               MPI_INT,		 //parentIndex
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //_packedRecords0
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               1,		 //hasToHoldDataForMasterWorkerCommunication
               DIMENSIONS_TIMES_TWO,		 //faceDataExchangeCounter
               1,		 //parentIndex
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //_packedRecords0
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._hasToHoldDataForMasterWorkerCommunication))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._faceDataExchangeCounter[0]))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._limiterStatus))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[32] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescriptionPacked::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescriptionPacked::Datatype );
         MPI_Type_free( &ADERDGCellDescriptionPacked::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "
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
               msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "
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
                  msg << "testing for finished send task for exahype::records::ADERDGCellDescriptionPacked "
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
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescriptionPacked",
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
      
      
      
      void exahype::records::ADERDGCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node "
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
               msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node "
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
                  msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: "
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
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::ADERDGCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      
      
   #endif
   
   
   
#elif !defined(Parallel)
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _riemannSolvePerformed(riemannSolvePerformed),
   _isInside(isInside),
   _parentIndex(parentIndex),
   _type(type),
   _refinementEvent(refinementEvent),
   _level(level),
   _offset(offset),
   _size(size),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _limiterStatus(limiterStatus),
   _mergedLimiterStatus(mergedLimiterStatus),
   _compressionState(compressionState),
   _bytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution),
   _bytesPerDoFInSolution(bytesPerDoFInSolution),
   _bytesPerDoFInUpdate(bytesPerDoFInUpdate),
   _bytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor),
   _bytesPerDoFInFluctuation(bytesPerDoFInFluctuation) {
      
   }
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription() {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords._riemannSolvePerformed, persistentRecords._isInside, persistentRecords._parentIndex, persistentRecords._type, persistentRecords._refinementEvent, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._limiterStatus, persistentRecords._mergedLimiterStatus, persistentRecords._compressionState, persistentRecords._bytesPerDoFInPreviousSolution, persistentRecords._bytesPerDoFInSolution, persistentRecords._bytesPerDoFInUpdate, persistentRecords._bytesPerDoFInExtrapolatedPredictor, persistentRecords._bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::ADERDGCellDescription(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, riemannSolvePerformed, isInside, parentIndex, type, refinementEvent, level, offset, size, correctorTimeStepSize, correctorTimeStamp, previousCorrectorTimeStepSize, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, limiterStatus, mergedLimiterStatus, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      
   }
   
   
   exahype::records::ADERDGCellDescription::~ADERDGCellDescription() { }
   
   std::string exahype::records::ADERDGCellDescription::toString(const CompressionState& param) {
      switch (param) {
         case Uncompressed: return "Uncompressed";
         case CurrentlyProcessed: return "CurrentlyProcessed";
         case Compressed: return "Compressed";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getCompressionStateMapping() {
      return "CompressionState(Uncompressed=0,CurrentlyProcessed=1,Compressed=2)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const LimiterStatus& param) {
      switch (param) {
         case Ok: return "Ok";
         case NeighbourIsNeighbourOfTroubledCell: return "NeighbourIsNeighbourOfTroubledCell";
         case NeighbourIsTroubledCell: return "NeighbourIsTroubledCell";
         case Troubled: return "Troubled";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getLimiterStatusMapping() {
      return "LimiterStatus(Ok=0,NeighbourIsNeighbourOfTroubledCell=1,NeighbourIsTroubledCell=2,Troubled=3)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const RefinementEvent& param) {
      switch (param) {
         case None: return "None";
         case ErasingChildrenRequested: return "ErasingChildrenRequested";
         case ErasingChildren: return "ErasingChildren";
         case ChangeChildrenToDescendantsRequested: return "ChangeChildrenToDescendantsRequested";
         case ChangeChildrenToDescendants: return "ChangeChildrenToDescendants";
         case RefiningRequested: return "RefiningRequested";
         case Refining: return "Refining";
         case DeaugmentingChildrenRequestedTriggered: return "DeaugmentingChildrenRequestedTriggered";
         case DeaugmentingChildrenRequested: return "DeaugmentingChildrenRequested";
         case DeaugmentingChildren: return "DeaugmentingChildren";
         case AugmentingRequested: return "AugmentingRequested";
         case Augmenting: return "Augmenting";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getRefinementEventMapping() {
      return "RefinementEvent(None=0,ErasingChildrenRequested=1,ErasingChildren=2,ChangeChildrenToDescendantsRequested=3,ChangeChildrenToDescendants=4,RefiningRequested=5,Refining=6,DeaugmentingChildrenRequestedTriggered=7,DeaugmentingChildrenRequested=8,DeaugmentingChildren=9,AugmentingRequested=10,Augmenting=11)";
   }
   std::string exahype::records::ADERDGCellDescription::toString(const Type& param) {
      switch (param) {
         case Erased: return "Erased";
         case Ancestor: return "Ancestor";
         case EmptyAncestor: return "EmptyAncestor";
         case Cell: return "Cell";
         case Descendant: return "Descendant";
         case EmptyDescendant: return "EmptyDescendant";
      }
      return "undefined";
   }
   
   std::string exahype::records::ADERDGCellDescription::getTypeMapping() {
      return "Type(Erased=0,Ancestor=1,EmptyAncestor=2,Cell=3,Descendant=4,EmptyDescendant=5)";
   }
   
   
   std::string exahype::records::ADERDGCellDescription::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescription::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "riemannSolvePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getRiemannSolvePerformed(i) << ",";
   }
   out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "isInside:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getIsInside(i) << ",";
   }
   out << getIsInside(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "limiterStatus:" << toString(getLimiterStatus());
      out << ",";
      out << "mergedLimiterStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getMergedLimiterStatus(i) << ",";
   }
   out << getMergedLimiterStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescription::PersistentRecords exahype::records::ADERDGCellDescription::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescriptionPacked exahype::records::ADERDGCellDescription::convert() const{
      return ADERDGCellDescriptionPacked(
         getSolverNumber(),
         getRiemannSolvePerformed(),
         getIsInside(),
         getParentIndex(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getLimiterStatus(),
         getMergedLimiterStatus(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescription::_log( "exahype::records::ADERDGCellDescription" );
      
      MPI_Datatype exahype::records::ADERDGCellDescription::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescription::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescription::initDatatype() {
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            const int Attributes = 40;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //riemannSolvePerformed
               MPI_INT,		 //isInside
               MPI_INT,		 //parentIndex
               MPI_INT,		 //type
               MPI_INT,		 //refinementEvent
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //compressionState
               MPI_INT,		 //bytesPerDoFInPreviousSolution
               MPI_INT,		 //bytesPerDoFInSolution
               MPI_INT,		 //bytesPerDoFInUpdate
               MPI_INT,		 //bytesPerDoFInExtrapolatedPredictor
               MPI_INT,		 //bytesPerDoFInFluctuation
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               DIMENSIONS_TIMES_TWO,		 //riemannSolvePerformed
               DIMENSIONS_TIMES_TWO,		 //isInside
               1,		 //parentIndex
               1,		 //type
               1,		 //refinementEvent
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //compressionState
               1,		 //bytesPerDoFInPreviousSolution
               1,		 //bytesPerDoFInSolution
               1,		 //bytesPerDoFInUpdate
               1,		 //bytesPerDoFInExtrapolatedPredictor
               1,		 //bytesPerDoFInFluctuation
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._riemannSolvePerformed))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._isInside))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._limiterStatus))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[32] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[33] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[34] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[35] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[36] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[37] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[38] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))), 		&disp[39] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::Datatype );
            MPI_Type_commit( &ADERDGCellDescription::Datatype );
            
         }
         {
            ADERDGCellDescription dummyADERDGCellDescription[2];
            
            const int Attributes = 40;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //riemannSolvePerformed
               MPI_INT,		 //isInside
               MPI_INT,		 //parentIndex
               MPI_INT,		 //type
               MPI_INT,		 //refinementEvent
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //compressionState
               MPI_INT,		 //bytesPerDoFInPreviousSolution
               MPI_INT,		 //bytesPerDoFInSolution
               MPI_INT,		 //bytesPerDoFInUpdate
               MPI_INT,		 //bytesPerDoFInExtrapolatedPredictor
               MPI_INT,		 //bytesPerDoFInFluctuation
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               DIMENSIONS_TIMES_TWO,		 //riemannSolvePerformed
               DIMENSIONS_TIMES_TWO,		 //isInside
               1,		 //parentIndex
               1,		 //type
               1,		 //refinementEvent
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //compressionState
               1,		 //bytesPerDoFInPreviousSolution
               1,		 //bytesPerDoFInSolution
               1,		 //bytesPerDoFInUpdate
               1,		 //bytesPerDoFInExtrapolatedPredictor
               1,		 //bytesPerDoFInFluctuation
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._riemannSolvePerformed))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._isInside))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._parentIndex))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._type))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._refinementEvent))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._level))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._offset[0]))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._size[0]))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStepSize))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._correctorTimeStamp))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStepSize))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._predictorTimeStamp))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solution))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionAverages))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionCompressed))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolution))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionAverages))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._previousSolutionCompressed))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._update))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateAverages))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._updateCompressed))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictor))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuation))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationAverages))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._fluctuationCompressed))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMin))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._solutionMax))), 		&disp[30] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._limiterStatus))), 		&disp[31] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[32] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._compressionState))), 		&disp[33] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInPreviousSolution))), 		&disp[34] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInSolution))), 		&disp[35] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInUpdate))), 		&disp[36] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInExtrapolatedPredictor))), 		&disp[37] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[0]._persistentRecords._bytesPerDoFInFluctuation))), 		&disp[38] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescription[1]._persistentRecords._solverNumber))), 		&disp[39] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescription::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescription::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescription::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescription::Datatype );
         MPI_Type_free( &ADERDGCellDescription::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::ADERDGCellDescription "
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
               msg << "was not able to send message exahype::records::ADERDGCellDescription "
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
                  msg << "testing for finished send task for exahype::records::ADERDGCellDescription "
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
                  "exahype::records::ADERDGCellDescription",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescription",
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
      
      
      
      void exahype::records::ADERDGCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::ADERDGCellDescription from node "
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
               msg << "failed to start to receive exahype::records::ADERDGCellDescription from node "
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
                  msg << "testing for finished receive task for exahype::records::ADERDGCellDescription failed: "
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
                  "exahype::records::ADERDGCellDescription",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescription",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::ADERDGCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      
      
   #endif
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords() {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _solverNumber(solverNumber),
   _parentIndex(parentIndex),
   _level(level),
   _offset(offset),
   _size(size),
   _correctorTimeStepSize(correctorTimeStepSize),
   _correctorTimeStamp(correctorTimeStamp),
   _previousCorrectorTimeStepSize(previousCorrectorTimeStepSize),
   _predictorTimeStepSize(predictorTimeStepSize),
   _predictorTimeStamp(predictorTimeStamp),
   _solution(solution),
   _solutionAverages(solutionAverages),
   _solutionCompressed(solutionCompressed),
   _previousSolution(previousSolution),
   _previousSolutionAverages(previousSolutionAverages),
   _previousSolutionCompressed(previousSolutionCompressed),
   _update(update),
   _updateAverages(updateAverages),
   _updateCompressed(updateCompressed),
   _extrapolatedPredictor(extrapolatedPredictor),
   _extrapolatedPredictorAverages(extrapolatedPredictorAverages),
   _extrapolatedPredictorCompressed(extrapolatedPredictorCompressed),
   _fluctuation(fluctuation),
   _fluctuationAverages(fluctuationAverages),
   _fluctuationCompressed(fluctuationCompressed),
   _solutionMin(solutionMin),
   _solutionMax(solutionMax),
   _limiterStatus(limiterStatus),
   _mergedLimiterStatus(mergedLimiterStatus) {
      setRiemannSolvePerformed(riemannSolvePerformed);
      setIsInside(isInside);
      setType(type);
      setRefinementEvent(refinementEvent);
      setCompressionState(compressionState);
      setBytesPerDoFInPreviousSolution(bytesPerDoFInPreviousSolution);
      setBytesPerDoFInSolution(bytesPerDoFInSolution);
      setBytesPerDoFInUpdate(bytesPerDoFInUpdate);
      setBytesPerDoFInExtrapolatedPredictor(bytesPerDoFInExtrapolatedPredictor);
      setBytesPerDoFInFluctuation(bytesPerDoFInFluctuation);
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked() {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._solverNumber, persistentRecords.getRiemannSolvePerformed(), persistentRecords.getIsInside(), persistentRecords._parentIndex, persistentRecords.getType(), persistentRecords.getRefinementEvent(), persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._correctorTimeStepSize, persistentRecords._correctorTimeStamp, persistentRecords._previousCorrectorTimeStepSize, persistentRecords._predictorTimeStepSize, persistentRecords._predictorTimeStamp, persistentRecords._solution, persistentRecords._solutionAverages, persistentRecords._solutionCompressed, persistentRecords._previousSolution, persistentRecords._previousSolutionAverages, persistentRecords._previousSolutionCompressed, persistentRecords._update, persistentRecords._updateAverages, persistentRecords._updateCompressed, persistentRecords._extrapolatedPredictor, persistentRecords._extrapolatedPredictorAverages, persistentRecords._extrapolatedPredictorCompressed, persistentRecords._fluctuation, persistentRecords._fluctuationAverages, persistentRecords._fluctuationCompressed, persistentRecords._solutionMin, persistentRecords._solutionMax, persistentRecords._limiterStatus, persistentRecords._mergedLimiterStatus, persistentRecords.getCompressionState(), persistentRecords.getBytesPerDoFInPreviousSolution(), persistentRecords.getBytesPerDoFInSolution(), persistentRecords.getBytesPerDoFInUpdate(), persistentRecords.getBytesPerDoFInExtrapolatedPredictor(), persistentRecords.getBytesPerDoFInFluctuation()) {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::ADERDGCellDescriptionPacked(const int& solverNumber, const std::bitset<DIMENSIONS_TIMES_TWO>& riemannSolvePerformed, const std::bitset<DIMENSIONS_TIMES_TWO>& isInside, const int& parentIndex, const Type& type, const RefinementEvent& refinementEvent, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const double& correctorTimeStepSize, const double& correctorTimeStamp, const double& previousCorrectorTimeStepSize, const double& predictorTimeStepSize, const double& predictorTimeStamp, const int& solution, const int& solutionAverages, const int& solutionCompressed, const int& previousSolution, const int& previousSolutionAverages, const int& previousSolutionCompressed, const int& update, const int& updateAverages, const int& updateCompressed, const int& extrapolatedPredictor, const int& extrapolatedPredictorAverages, const int& extrapolatedPredictorCompressed, const int& fluctuation, const int& fluctuationAverages, const int& fluctuationCompressed, const int& solutionMin, const int& solutionMax, const LimiterStatus& limiterStatus, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,LimiterStatus>& mergedLimiterStatus, const CompressionState& compressionState, const int& bytesPerDoFInPreviousSolution, const int& bytesPerDoFInSolution, const int& bytesPerDoFInUpdate, const int& bytesPerDoFInExtrapolatedPredictor, const int& bytesPerDoFInFluctuation):
   _persistentRecords(solverNumber, riemannSolvePerformed, isInside, parentIndex, type, refinementEvent, level, offset, size, correctorTimeStepSize, correctorTimeStamp, previousCorrectorTimeStepSize, predictorTimeStepSize, predictorTimeStamp, solution, solutionAverages, solutionCompressed, previousSolution, previousSolutionAverages, previousSolutionCompressed, update, updateAverages, updateCompressed, extrapolatedPredictor, extrapolatedPredictorAverages, extrapolatedPredictorCompressed, fluctuation, fluctuationAverages, fluctuationCompressed, solutionMin, solutionMax, limiterStatus, mergedLimiterStatus, compressionState, bytesPerDoFInPreviousSolution, bytesPerDoFInSolution, bytesPerDoFInUpdate, bytesPerDoFInExtrapolatedPredictor, bytesPerDoFInFluctuation) {
      if ((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 >= (8 * sizeof(int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((DIMENSIONS_TIMES_TWO+DIMENSIONS_TIMES_TWO+24 < (8 * sizeof(int))));
      
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::~ADERDGCellDescriptionPacked() { }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const Type& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getTypeMapping() {
      return exahype::records::ADERDGCellDescription::getTypeMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const RefinementEvent& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getRefinementEventMapping() {
      return exahype::records::ADERDGCellDescription::getRefinementEventMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const LimiterStatus& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getLimiterStatusMapping() {
      return exahype::records::ADERDGCellDescription::getLimiterStatusMapping();
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString(const CompressionState& param) {
      return exahype::records::ADERDGCellDescription::toString(param);
   }
   
   std::string exahype::records::ADERDGCellDescriptionPacked::getCompressionStateMapping() {
      return exahype::records::ADERDGCellDescription::getCompressionStateMapping();
   }
   
   
   
   std::string exahype::records::ADERDGCellDescriptionPacked::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::ADERDGCellDescriptionPacked::toString (std::ostream& out) const {
      out << "("; 
      out << "solverNumber:" << getSolverNumber();
      out << ",";
      out << "riemannSolvePerformed:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getRiemannSolvePerformed(i) << ",";
   }
   out << getRiemannSolvePerformed(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "isInside:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getIsInside(i) << ",";
   }
   out << getIsInside(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "parentIndex:" << getParentIndex();
      out << ",";
      out << "type:" << toString(getType());
      out << ",";
      out << "refinementEvent:" << toString(getRefinementEvent());
      out << ",";
      out << "level:" << getLevel();
      out << ",";
      out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
      out << ",";
      out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
      out << ",";
      out << "correctorTimeStepSize:" << getCorrectorTimeStepSize();
      out << ",";
      out << "correctorTimeStamp:" << getCorrectorTimeStamp();
      out << ",";
      out << "previousCorrectorTimeStepSize:" << getPreviousCorrectorTimeStepSize();
      out << ",";
      out << "predictorTimeStepSize:" << getPredictorTimeStepSize();
      out << ",";
      out << "predictorTimeStamp:" << getPredictorTimeStamp();
      out << ",";
      out << "solution:" << getSolution();
      out << ",";
      out << "solutionAverages:" << getSolutionAverages();
      out << ",";
      out << "solutionCompressed:" << getSolutionCompressed();
      out << ",";
      out << "previousSolution:" << getPreviousSolution();
      out << ",";
      out << "previousSolutionAverages:" << getPreviousSolutionAverages();
      out << ",";
      out << "previousSolutionCompressed:" << getPreviousSolutionCompressed();
      out << ",";
      out << "update:" << getUpdate();
      out << ",";
      out << "updateAverages:" << getUpdateAverages();
      out << ",";
      out << "updateCompressed:" << getUpdateCompressed();
      out << ",";
      out << "extrapolatedPredictor:" << getExtrapolatedPredictor();
      out << ",";
      out << "extrapolatedPredictorAverages:" << getExtrapolatedPredictorAverages();
      out << ",";
      out << "extrapolatedPredictorCompressed:" << getExtrapolatedPredictorCompressed();
      out << ",";
      out << "fluctuation:" << getFluctuation();
      out << ",";
      out << "fluctuationAverages:" << getFluctuationAverages();
      out << ",";
      out << "fluctuationCompressed:" << getFluctuationCompressed();
      out << ",";
      out << "solutionMin:" << getSolutionMin();
      out << ",";
      out << "solutionMax:" << getSolutionMax();
      out << ",";
      out << "limiterStatus:" << toString(getLimiterStatus());
      out << ",";
      out << "mergedLimiterStatus:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getMergedLimiterStatus(i) << ",";
   }
   out << getMergedLimiterStatus(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "compressionState:" << toString(getCompressionState());
      out << ",";
      out << "bytesPerDoFInPreviousSolution:" << getBytesPerDoFInPreviousSolution();
      out << ",";
      out << "bytesPerDoFInSolution:" << getBytesPerDoFInSolution();
      out << ",";
      out << "bytesPerDoFInUpdate:" << getBytesPerDoFInUpdate();
      out << ",";
      out << "bytesPerDoFInExtrapolatedPredictor:" << getBytesPerDoFInExtrapolatedPredictor();
      out << ",";
      out << "bytesPerDoFInFluctuation:" << getBytesPerDoFInFluctuation();
      out <<  ")";
   }
   
   
   exahype::records::ADERDGCellDescriptionPacked::PersistentRecords exahype::records::ADERDGCellDescriptionPacked::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::ADERDGCellDescription exahype::records::ADERDGCellDescriptionPacked::convert() const{
      return ADERDGCellDescription(
         getSolverNumber(),
         getRiemannSolvePerformed(),
         getIsInside(),
         getParentIndex(),
         getType(),
         getRefinementEvent(),
         getLevel(),
         getOffset(),
         getSize(),
         getCorrectorTimeStepSize(),
         getCorrectorTimeStamp(),
         getPreviousCorrectorTimeStepSize(),
         getPredictorTimeStepSize(),
         getPredictorTimeStamp(),
         getSolution(),
         getSolutionAverages(),
         getSolutionCompressed(),
         getPreviousSolution(),
         getPreviousSolutionAverages(),
         getPreviousSolutionCompressed(),
         getUpdate(),
         getUpdateAverages(),
         getUpdateCompressed(),
         getExtrapolatedPredictor(),
         getExtrapolatedPredictorAverages(),
         getExtrapolatedPredictorCompressed(),
         getFluctuation(),
         getFluctuationAverages(),
         getFluctuationCompressed(),
         getSolutionMin(),
         getSolutionMax(),
         getLimiterStatus(),
         getMergedLimiterStatus(),
         getCompressionState(),
         getBytesPerDoFInPreviousSolution(),
         getBytesPerDoFInSolution(),
         getBytesPerDoFInUpdate(),
         getBytesPerDoFInExtrapolatedPredictor(),
         getBytesPerDoFInFluctuation()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::ADERDGCellDescriptionPacked::_log( "exahype::records::ADERDGCellDescriptionPacked" );
      
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::Datatype = 0;
      MPI_Datatype exahype::records::ADERDGCellDescriptionPacked::FullDatatype = 0;
      
      
      void exahype::records::ADERDGCellDescriptionPacked::initDatatype() {
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            const int Attributes = 31;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //parentIndex
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //_packedRecords0
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               1,		 //parentIndex
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //_packedRecords0
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._limiterStatus))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[30] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::Datatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::Datatype );
            
         }
         {
            ADERDGCellDescriptionPacked dummyADERDGCellDescriptionPacked[2];
            
            const int Attributes = 31;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //solverNumber
               MPI_INT,		 //parentIndex
               MPI_INT,		 //level
               MPI_DOUBLE,		 //offset
               MPI_DOUBLE,		 //size
               MPI_DOUBLE,		 //correctorTimeStepSize
               MPI_DOUBLE,		 //correctorTimeStamp
               MPI_DOUBLE,		 //previousCorrectorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStepSize
               MPI_DOUBLE,		 //predictorTimeStamp
               MPI_INT,		 //solution
               MPI_INT,		 //solutionAverages
               MPI_INT,		 //solutionCompressed
               MPI_INT,		 //previousSolution
               MPI_INT,		 //previousSolutionAverages
               MPI_INT,		 //previousSolutionCompressed
               MPI_INT,		 //update
               MPI_INT,		 //updateAverages
               MPI_INT,		 //updateCompressed
               MPI_INT,		 //extrapolatedPredictor
               MPI_INT,		 //extrapolatedPredictorAverages
               MPI_INT,		 //extrapolatedPredictorCompressed
               MPI_INT,		 //fluctuation
               MPI_INT,		 //fluctuationAverages
               MPI_INT,		 //fluctuationCompressed
               MPI_INT,		 //solutionMin
               MPI_INT,		 //solutionMax
               MPI_INT,		 //limiterStatus
               MPI_INT,		 //mergedLimiterStatus
               MPI_INT,		 //_packedRecords0
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //solverNumber
               1,		 //parentIndex
               1,		 //level
               DIMENSIONS,		 //offset
               DIMENSIONS,		 //size
               1,		 //correctorTimeStepSize
               1,		 //correctorTimeStamp
               1,		 //previousCorrectorTimeStepSize
               1,		 //predictorTimeStepSize
               1,		 //predictorTimeStamp
               1,		 //solution
               1,		 //solutionAverages
               1,		 //solutionCompressed
               1,		 //previousSolution
               1,		 //previousSolutionAverages
               1,		 //previousSolutionCompressed
               1,		 //update
               1,		 //updateAverages
               1,		 //updateCompressed
               1,		 //extrapolatedPredictor
               1,		 //extrapolatedPredictorAverages
               1,		 //extrapolatedPredictorCompressed
               1,		 //fluctuation
               1,		 //fluctuationAverages
               1,		 //fluctuationCompressed
               1,		 //solutionMin
               1,		 //solutionMax
               1,		 //limiterStatus
               DIMENSIONS_TIMES_TWO,		 //mergedLimiterStatus
               1,		 //_packedRecords0
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._parentIndex))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStepSize))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._correctorTimeStamp))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousCorrectorTimeStepSize))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStepSize))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._predictorTimeStamp))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionAverages))), 		&disp[11] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionCompressed))), 		&disp[12] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolution))), 		&disp[13] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionAverages))), 		&disp[14] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._previousSolutionCompressed))), 		&disp[15] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._update))), 		&disp[16] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateAverages))), 		&disp[17] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._updateCompressed))), 		&disp[18] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor))), 		&disp[19] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorAverages))), 		&disp[20] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictorCompressed))), 		&disp[21] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuation))), 		&disp[22] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationAverages))), 		&disp[23] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._fluctuationCompressed))), 		&disp[24] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMin))), 		&disp[25] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._solutionMax))), 		&disp[26] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._limiterStatus))), 		&disp[27] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._mergedLimiterStatus[0]))), 		&disp[28] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[0]._persistentRecords._packedRecords0))), 		&disp[29] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyADERDGCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[30] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ADERDGCellDescriptionPacked::FullDatatype );
            MPI_Type_commit( &ADERDGCellDescriptionPacked::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::ADERDGCellDescriptionPacked::shutdownDatatype() {
         MPI_Type_free( &ADERDGCellDescriptionPacked::Datatype );
         MPI_Type_free( &ADERDGCellDescriptionPacked::FullDatatype );
         
      }
      
      void exahype::records::ADERDGCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "
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
               msg << "was not able to send message exahype::records::ADERDGCellDescriptionPacked "
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
                  msg << "testing for finished send task for exahype::records::ADERDGCellDescriptionPacked "
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
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescriptionPacked",
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
      
      
      
      void exahype::records::ADERDGCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node "
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
               msg << "failed to start to receive exahype::records::ADERDGCellDescriptionPacked from node "
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
                  msg << "testing for finished receive task for exahype::records::ADERDGCellDescriptionPacked failed: "
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
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::ADERDGCellDescriptionPacked",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::ADERDGCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
      
      
   #endif
   
   
   

#endif


