#include "SHMLockTask.h"
#include "SHMController.h"


#include <iostream>
#include <assert.h>


shminvade::SHMLockTask::SHMLockTask(pid_t pid_t):
  _pid_t( pid_t ) {
  // TBB complains and clarifies that affinities are ignored for enqueued
  // tasks
  // set_affinity( _pid_t );
}


void shminvade::SHMLockTask::reenqueue() {
  #if SHM_INVADE_DEBUG>=8
  std::cout << SHM_DEBUG_PREFIX <<  "Lock task for thread id " << _pid_t <<
    " has been triggered to reenqueue on thread " << syscall (__NR_gettid) << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  if (SHMController::getInstance().getThreadTableEntry(_pid_t).type == SHMController::ThreadType::Shutdown) {
    #if SHM_INVADE_DEBUG>=8
    std::cout << SHM_DEBUG_PREFIX <<  "Controller seems to be down already, so stop lock task too (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }
  else {
    tbb::task &t = *new(tbb::task::allocate_root()) SHMLockTask(_pid_t);
    tbb::task::enqueue(t);
  }
}


void shminvade::SHMLockTask::terminate() {
  SHMController::ThreadTable::accessor a;
  SHMController::getInstance()._threads.find(a,_pid_t);
  assert( a->second.numberOfExistingLockTasks>0 );
  a->second.numberOfExistingLockTasks--;
}



tbb::task* shminvade::SHMLockTask::execute() {
  const pid_t currentThreadId = (pid_t) syscall (__NR_gettid);

  if ( SHMController::getInstance()._threads.count(currentThreadId)==0 ) {
    SHMController::getInstance().registerNewThread(currentThreadId);
  }

  if ( currentThreadId!=_pid_t ) {
    reenqueue();
    return nullptr;
  }
  else {
    SHMController::ThreadType state = SHMController::getInstance().getThreadTableEntry(_pid_t).type;
    switch (state) {
      case SHMController::ThreadType::Master:
        // there should be no lock tasks on the master so let this one die
        #if SHM_INVADE_DEBUG>=8
        std::cout << SHM_DEBUG_PREFIX <<  "thread " << _pid_t << " is the master so remove lock thread (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
        #endif
        terminate();
        return nullptr;
      case SHMController::ThreadType::ExclusivelyOwned:
        // we own it so let this one die
        #if SHM_INVADE_DEBUG>=8
        std::cout << SHM_DEBUG_PREFIX <<  "thread " << _pid_t << " is exclusively owned so remove lock thread (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
        #endif
        terminate();
        return nullptr;
      case SHMController::ThreadType::NotOwned:
        #if SHM_INVADE_DEBUG>=8
        std::cout << SHM_DEBUG_PREFIX <<  "thread " << _pid_t << " should not be used so yield immediately and make lock task sleep (line:" << __LINE__ << ",file: " << __FILE__ << ")" << std::endl;
        #endif
        __TBB_Yield();
        #ifdef SHM_SLEEP
        sleep(SHM_SLEEP);
        #endif
        reenqueue();
        return nullptr;
      case SHMController::ThreadType::Shutdown:
        #if SHM_INVADE_DEBUG>=8
        std::cout << SHM_DEBUG_PREFIX <<  "thread " << _pid_t << " is marked to shut down so remove lock thread (line:" << __LINE__  << ",file: " << __FILE__ << ")" << std::endl;
        #endif
        // we should die anyway
        terminate();
        return nullptr;
    }
  }
  return nullptr;
}
