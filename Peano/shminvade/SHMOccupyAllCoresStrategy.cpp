#include "SHMOccupyAllCoresStrategy.h"
#include "SHMController.h"
#include "SHMMacros.h"


#include <iostream>


shminvade::SHMOccupyAllCoresStrategy::SHMOccupyAllCoresStrategy() {
}


shminvade::SHMOccupyAllCoresStrategy::~SHMOccupyAllCoresStrategy() {
}


std::set<pid_t> shminvade::SHMOccupyAllCoresStrategy::invadeThreads(int wantedNumberOfThreads) {
  std::set<pid_t> bookedCores;

  if (SHMController::getInstance()._switchedOn) {
    typedef tbb::spin_mutex CoreTraversalMutex;
    static CoreTraversalMutex invadeMutex;
    CoreTraversalMutex::scoped_lock lock(invadeMutex);

    for (auto p: SHMController::getInstance()._threads) {
      SHMController::ThreadTable::accessor a;
      SHMController::getInstance()._threads.find(a,p.first);
      if (
        wantedNumberOfThreads>0
        &&
        ( a->second.type==SHMController::ThreadType::NotOwned )
      ) {
        wantedNumberOfThreads--;
        a->second.type = SHMController::ThreadType::ExclusivelyOwned;
        bookedCores.insert(a->first);
        #if SHM_INVADE_DEBUG>=4
        std::cout << SHM_DEBUG_PREFIX <<  "invade core " << a->first << " (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
        #endif
      }
    }

    lock.release();
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "invaded " << bookedCores.size() << " core(s) in total with " << wantedNumberOfThreads << " open requests (line:" << __LINE__  << ",file: " << __FILE__ <<  ")" << std::endl;
    #endif
  }

  return bookedCores;
}


void shminvade::SHMOccupyAllCoresStrategy::cleanUp() {
}


void shminvade::SHMOccupyAllCoresStrategy::retreat(const std::set<pid_t>& threadIds) {
  for (auto p: threadIds) {
    SHMController::getInstance().retreatFromCore(p);
  }
}
