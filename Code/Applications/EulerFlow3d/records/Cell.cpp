#include "EulerFlow3d/records/Cell.h"

#if !defined(Debug) && !defined(Parallel) && defined(SharedMemoryParallelisation)
   exahype::records::Cell::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
   _cellDescriptionsIndex(cellDescriptionsIndex),
   _isInside(isInside),
   _state(state),
   _evenFlags(evenFlags),
   _accessNumber(accessNumber),
   _numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
   _numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
      
   }
   
   exahype::records::Cell::Cell() {
      
   }
   
   
   exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
      
   }
   
   
   exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
   _persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
      
   }
   
   
   exahype::records::Cell::~Cell() { }
   
   std::string exahype::records::Cell::toString(const State& param) {
      switch (param) {
         case Leaf: return "Leaf";
         case Refined: return "Refined";
         case Root: return "Root";
      }
      return "undefined";
   }
   
   std::string exahype::records::Cell::getStateMapping() {
      return "State(Leaf=0,Refined=1,Root=2)";
   }
   
   
   std::string exahype::records::Cell::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::Cell::toString (std::ostream& out) const {
      out << "("; 
      out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
      out << ",";
      out << "isInside:" << getIsInside();
      out << ",";
      out << "state:" << toString(getState());
      out << ",";
      out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
      out << ",";
      out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
      out << ",";
      out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
      out <<  ")";
   }
   
   
   exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::CellPacked exahype::records::Cell::convert() const{
      return CellPacked(
         getCellDescriptionsIndex(),
         getIsInside(),
         getState(),
         getEvenFlags(),
         getAccessNumber(),
         getNumberOfLoadsFromInputStream(),
         getNumberOfStoresToOutputStream()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );
      
      MPI_Datatype exahype::records::Cell::Datatype = 0;
      MPI_Datatype exahype::records::Cell::FullDatatype = 0;
      
      
      void exahype::records::Cell::initDatatype() {
         {
            Cell dummyCell[2];
            
            const int Attributes = 4;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //cellDescriptionsIndex
               MPI_CHAR,		 //isInside
               MPI_INT,		 //state
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //cellDescriptionsIndex
               1,		 //isInside
               1,		 //state
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[3] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
            MPI_Type_commit( &Cell::Datatype );
            
         }
         {
            Cell dummyCell[2];
            
            const int Attributes = 8;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //cellDescriptionsIndex
               MPI_CHAR,		 //isInside
               MPI_INT,		 //state
               MPI_INT,		 //evenFlags
               MPI_SHORT,		 //accessNumber
               MPI_INT,		 //numberOfLoadsFromInputStream
               MPI_INT,		 //numberOfStoresToOutputStream
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //cellDescriptionsIndex
               1,		 //isInside
               1,		 //state
               DIMENSIONS,		 //evenFlags
               DIMENSIONS_TIMES_TWO,		 //accessNumber
               1,		 //numberOfLoadsFromInputStream
               1,		 //numberOfStoresToOutputStream
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[7] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
            MPI_Type_commit( &Cell::FullDatatype );
            
         }
         
      }
      
      
      void exahype::records::Cell::shutdownDatatype() {
         MPI_Type_free( &Cell::Datatype );
         MPI_Type_free( &Cell::FullDatatype );
         
      }
      
      void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         _senderDestinationRank = destination;
         
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::Cell "
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
            msg << "was not able to send message exahype::records::Cell "
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
               msg << "testing for finished send task for exahype::records::Cell "
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
               "exahype::records::Cell",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::Cell",
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
   
   
   
   void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         _senderDestinationRank = status.MPI_SOURCE;
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::Cell from node "
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
            msg << "failed to start to receive exahype::records::Cell from node "
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
               msg << "testing for finished receive task for exahype::records::Cell failed: "
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
               "exahype::records::Cell",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::Cell",
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
   
   
   
   bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int exahype::records::Cell::getSenderRank() const {
      assertion( _senderDestinationRank!=-1 );
      return _senderDestinationRank;
      
   }
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
   if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((DIMENSIONS+3 < (8 * sizeof(int))));
   
}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
   setIsInside(isInside);
   setState(state);
   setEvenFlags(evenFlags);
   if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((DIMENSIONS+3 < (8 * sizeof(int))));
   
}

exahype::records::CellPacked::CellPacked() {
   if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((DIMENSIONS+3 < (8 * sizeof(int))));
   
}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
   if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((DIMENSIONS+3 < (8 * sizeof(int))));
   
}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
   if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
      std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
      std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
      std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
   }
   assertion((DIMENSIONS+3 < (8 * sizeof(int))));
   
}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
   return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
   return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
   out << ",";
   out << "isInside:" << getIsInside();
   out << ",";
   out << "state:" << toString(getState());
   out << ",";
   out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
   out << ",";
   out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
   out << ",";
   out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
   out << ",";
   out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
   out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
   return Cell(
      getCellDescriptionsIndex(),
      getIsInside(),
      getState(),
      getEvenFlags(),
      getAccessNumber(),
      getNumberOfLoadsFromInputStream(),
      getNumberOfStoresToOutputStream()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );
   
   MPI_Datatype exahype::records::CellPacked::Datatype = 0;
   MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;
   
   
   void exahype::records::CellPacked::initDatatype() {
      {
         CellPacked dummyCellPacked[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //cellDescriptionsIndex
            MPI_INT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //cellDescriptionsIndex
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
         MPI_Type_commit( &CellPacked::Datatype );
         
      }
      {
         CellPacked dummyCellPacked[2];
         
         const int Attributes = 6;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //cellDescriptionsIndex
            MPI_SHORT,		 //accessNumber
            MPI_INT,		 //numberOfLoadsFromInputStream
            MPI_INT,		 //numberOfStoresToOutputStream
            MPI_INT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //cellDescriptionsIndex
            DIMENSIONS_TIMES_TWO,		 //accessNumber
            1,		 //numberOfLoadsFromInputStream
            1,		 //numberOfStoresToOutputStream
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[5] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
         MPI_Type_commit( &CellPacked::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::CellPacked::shutdownDatatype() {
      MPI_Type_free( &CellPacked::Datatype );
      MPI_Type_free( &CellPacked::FullDatatype );
      
   }
   
   void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      _senderDestinationRank = destination;
      
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::CellPacked "
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
         msg << "was not able to send message exahype::records::CellPacked "
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
            msg << "testing for finished send task for exahype::records::CellPacked "
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
            "exahype::records::CellPacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   if (communicateSleep<0) {
   
      MPI_Status  status;
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
      _senderDestinationRank = status.MPI_SOURCE;
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive exahype::records::CellPacked from node "
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
         msg << "failed to start to receive exahype::records::CellPacked from node "
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
            msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
            "exahype::records::CellPacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
   assertion( _senderDestinationRank!=-1 );
   return _senderDestinationRank;
   
}
#endif



#elif !defined(Parallel) && defined(Debug) && !defined(SharedMemoryParallelisation)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
   case Leaf: return "Leaf";
   case Refined: return "Refined";
   case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
   getCellDescriptionsIndex(),
   getIsInside(),
   getState(),
   getLevel(),
   getEvenFlags(),
   getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
   {
      Cell dummyCell[2];
      
      const int Attributes = 5;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //cellDescriptionsIndex
         MPI_CHAR,		 //isInside
         MPI_INT,		 //state
         MPI_INT,		 //level
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //cellDescriptionsIndex
         1,		 //isInside
         1,		 //state
         1,		 //level
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[4] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
      MPI_Type_commit( &Cell::Datatype );
      
   }
   {
      Cell dummyCell[2];
      
      const int Attributes = 7;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //cellDescriptionsIndex
         MPI_CHAR,		 //isInside
         MPI_INT,		 //state
         MPI_INT,		 //level
         MPI_INT,		 //evenFlags
         MPI_SHORT,		 //accessNumber
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //cellDescriptionsIndex
         1,		 //isInside
         1,		 //state
         1,		 //level
         DIMENSIONS,		 //evenFlags
         DIMENSIONS_TIMES_TWO,		 //accessNumber
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[4] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[6] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
      MPI_Type_commit( &Cell::FullDatatype );
      
   }
   
}


void exahype::records::Cell::shutdownDatatype() {
   MPI_Type_free( &Cell::Datatype );
   MPI_Type_free( &Cell::FullDatatype );
   
}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
   _senderDestinationRank = destination;
   
   if (communicateSleep<0) {
   
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message exahype::records::Cell "
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
      msg << "was not able to send message exahype::records::Cell "
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
         msg << "testing for finished send task for exahype::records::Cell "
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
         "exahype::records::Cell",
         "send(int)", destination,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

   MPI_Status  status;
   const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
   _senderDestinationRank = status.MPI_SOURCE;
   if ( result != MPI_SUCCESS ) {
      std::ostringstream msg;
      msg << "failed to start to receive exahype::records::Cell from node "
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
      msg << "failed to start to receive exahype::records::Cell from node "
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
         msg << "testing for finished receive task for exahype::records::Cell failed: "
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
         "exahype::records::Cell",
         "receive(int)", source,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellDescriptionsIndex(cellDescriptionsIndex),
_level(level),
_accessNumber(accessNumber) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
   CellPacked dummyCellPacked[2];
   
   const int Attributes = 4;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //cellDescriptionsIndex
      MPI_INT,		 //level
      MPI_INT,		 //_packedRecords0
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //cellDescriptionsIndex
      1,		 //level
      1,		 //_packedRecords0
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[3] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
   MPI_Type_commit( &CellPacked::Datatype );
   
}
{
   CellPacked dummyCellPacked[2];
   
   const int Attributes = 5;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //cellDescriptionsIndex
      MPI_INT,		 //level
      MPI_SHORT,		 //accessNumber
      MPI_INT,		 //_packedRecords0
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //cellDescriptionsIndex
      1,		 //level
      DIMENSIONS_TIMES_TWO,		 //accessNumber
      1,		 //_packedRecords0
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[3] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[4] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
   MPI_Type_commit( &CellPacked::FullDatatype );
   
}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

   const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
   if  (result!=MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message exahype::records::CellPacked "
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
   msg << "was not able to send message exahype::records::CellPacked "
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
      msg << "testing for finished send task for exahype::records::CellPacked "
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
      "exahype::records::CellPacked",
      "send(int)", destination,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
   std::ostringstream msg;
   msg << "failed to start to receive exahype::records::CellPacked from node "
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
   msg << "failed to start to receive exahype::records::CellPacked from node "
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
      msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
      "exahype::records::CellPacked",
      "receive(int)", source,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && !defined(Debug) && !defined(SharedMemoryParallelisation)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 14;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[13] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellDescriptionsIndex(cellDescriptionsIndex),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate()) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 9;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[8] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif !defined(Debug) && !defined(Parallel) && !defined(SharedMemoryParallelisation)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 4;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[3] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 6;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[5] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellDescriptionsIndex(cellDescriptionsIndex),
_accessNumber(accessNumber) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 3;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[2] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 4;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_SHORT,		 //accessNumber
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[3] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && defined(SharedMemoryParallelisation) && defined(Debug)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 17;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[16] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_level(level),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate(), persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 14;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[13] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && defined(Debug) && !defined(SharedMemoryParallelisation)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 15;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[14] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellDescriptionsIndex(cellDescriptionsIndex),
_level(level),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate()) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 12;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[11] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && !defined(Debug) && defined(SharedMemoryParallelisation)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 16;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[15] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate(), persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
if ((DIMENSIONS+4 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+4 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 9;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[8] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 13;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[12] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif !defined(Parallel) && defined(SharedMemoryParallelisation) && defined(Debug)
exahype::records::Cell::PersistentRecords::PersistentRecords() {

}


exahype::records::Cell::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

exahype::records::Cell::Cell() {

}


exahype::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


exahype::records::Cell::Cell(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


exahype::records::Cell::~Cell() { }

std::string exahype::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string exahype::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string exahype::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::Cell::PersistentRecords exahype::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::CellPacked exahype::records::Cell::convert() const{
return CellPacked(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::Cell::_log( "exahype::records::Cell" );

MPI_Datatype exahype::records::Cell::Datatype = 0;
MPI_Datatype exahype::records::Cell::FullDatatype = 0;


void exahype::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 5;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[4] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 9;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[8] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void exahype::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void exahype::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::Cell "
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
msg << "was not able to send message exahype::records::Cell "
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
msg << "testing for finished send task for exahype::records::Cell "
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
"exahype::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



void exahype::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "failed to start to receive exahype::records::Cell from node "
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
msg << "testing for finished receive task for exahype::records::Cell failed: "
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
"exahype::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::Cell",
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



bool exahype::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


exahype::records::CellPacked::PersistentRecords::PersistentRecords() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellDescriptionsIndex(cellDescriptionsIndex),
_level(level),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}

exahype::records::CellPacked::CellPacked() {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellDescriptionsIndex, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::CellPacked(const int& cellDescriptionsIndex, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellDescriptionsIndex, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
if ((DIMENSIONS+3 >= (8 * sizeof(int)))) {
std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
std::cerr << "  Packed-Type: int hint-size no-of-bits;  " << std::endl << std::endl;
std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
}
assertion((DIMENSIONS+3 < (8 * sizeof(int))));

}


exahype::records::CellPacked::~CellPacked() { }

std::string exahype::records::CellPacked::toString(const State& param) {
return exahype::records::Cell::toString(param);
}

std::string exahype::records::CellPacked::getStateMapping() {
return exahype::records::Cell::getStateMapping();
}



std::string exahype::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void exahype::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellDescriptionsIndex:" << getCellDescriptionsIndex();
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


exahype::records::CellPacked::PersistentRecords exahype::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

exahype::records::Cell exahype::records::CellPacked::convert() const{
return Cell(
getCellDescriptionsIndex(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log exahype::records::CellPacked::_log( "exahype::records::CellPacked" );

MPI_Datatype exahype::records::CellPacked::Datatype = 0;
MPI_Datatype exahype::records::CellPacked::FullDatatype = 0;


void exahype::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 4;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[3] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 7;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellDescriptionsIndex
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //_packedRecords0
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellDescriptionsIndex
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellDescriptionsIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellDescriptionsIndex))), 		&disp[6] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void exahype::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void exahype::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
_senderDestinationRank = destination;

if (communicateSleep<0) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "was not able to send message exahype::records::CellPacked "
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
msg << "testing for finished send task for exahype::records::CellPacked "
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
"exahype::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



void exahype::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
if (communicateSleep<0) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "failed to start to receive exahype::records::CellPacked from node "
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
msg << "testing for finished receive task for exahype::records::CellPacked failed: "
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
"exahype::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"exahype::records::CellPacked",
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



bool exahype::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int exahype::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#endif


