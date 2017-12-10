#include "SHMController.h"
#include "SHMLockTask.h"


#include <iostream>
#include <thread>
#include <sstream>


shminvade::SHMController::SHMController():
  _switchedOn( true ) {
  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Create SHMController, maxCores=" << getMaxAvailableCores() <<
    " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  registerNewThread(
    (pid_t) syscall (__NR_gettid),
    ThreadType::Master
  );
}


shminvade::SHMController::~SHMController() {
  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "Destroy SHMController (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  shutdown();
}


shminvade::SHMController&  shminvade::SHMController::getInstance() {
  static SHMController singleton;
  return singleton;
}


void shminvade::SHMController::switchOn() {
  _switchedOn = true;
}


void shminvade::SHMController::switchOff() {
  //_switchedOn = false;
}


int shminvade::SHMController::getMaxAvailableCores() const {
  return std::thread::hardware_concurrency();
}


int shminvade::SHMController::getFreeCores() const {
  return getMaxAvailableCores() - getBookedCores();
}


int shminvade::SHMController::getBookedCores() const {
  int result = 1;

  for (auto p: _threads) {
    if ( getThreadTableEntry(p.first).type==ThreadType::ExclusivelyOwned) {
      result++;
    }
  }

  return result;
}


shminvade::SHMController::ThreadState shminvade::SHMController::getThreadTableEntry( pid_t pid ) const {
  // switch to non-const to ensure that noone else is modifying the entry at the very moment
  ThreadTable::const_accessor a;
  _threads.find(a,pid);
  return a->second;
}


void shminvade::SHMController::shutdown() {
  _switchedOn = false;

  #if SHM_INVADE_DEBUG>=2
  std::cout << SHM_DEBUG_PREFIX <<  "start to instruct all threads to shut down (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  for (auto p: _threads) {
    ThreadTable::accessor a;
    _threads.find(a,p.first);
    a->second.type = ThreadType::Shutdown;
  }

  // @todo We should check whether there are lock tasks remainidng. If not, cancel the sleep

  #if SHM_INVADE_DEBUG>=2
  std::cout << SHM_DEBUG_PREFIX <<  "wait for all lock threads to terminate (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  __TBB_Yield();
  #ifdef SHM_SLEEP
  sleep(SHM_SLEEP);
  //sleep(getMaxAvailableCores() * SHM_SLEEP);
  #endif

  #if SHM_INVADE_DEBUG>=2
  std::cout << SHM_DEBUG_PREFIX <<  "assume all threads have terminated (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif
}


void shminvade::SHMController::retreatFromCore( pid_t pid ) {
  ThreadTable::accessor a;
  _threads.find(a,pid);
  if ( a->second.type==ThreadType::ExclusivelyOwned ) {
    a->second.type = ThreadType::NotOwned;
  }

  if (
    a->second.type!=ThreadType::Shutdown
    &&
    a->second.numberOfExistingLockTasks==0
  ) {
    a->second.numberOfExistingLockTasks++;
    tbb::task &t = *new(tbb::task::allocate_root()) SHMLockTask(pid);
    tbb::task::enqueue(t);
    #if SHM_INVADE_DEBUG>=4
    std::cout << SHM_DEBUG_PREFIX <<  "retreat from core " << pid << " and issue new lock task (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
    #endif
  }
}


void shminvade::SHMController::retreatFromAllCores() {
  for (auto p: _threads) {
    if (getThreadTableEntry(p.first).type==ThreadType::ExclusivelyOwned) {
      retreatFromCore(p.first);
    }
  }
}


void shminvade::SHMController::registerNewThread(int threadId, ThreadType initialType) {
  ThreadState newThread(initialType);
  _threads.insert( std::pair<pid_t, ThreadState>(threadId,newThread) );

  #if SHM_INVADE_DEBUG>=1
  std::cout << SHM_DEBUG_PREFIX <<  "register new thread " << threadId << " as " << toString(newThread) << " (line:" << __LINE__ << ",file:" << __FILE__ << ")" << std::endl;
  #endif

  retreatFromCore(threadId);
}


std::string shminvade::SHMController::toString( ThreadState state ) {
  std::ostringstream msg;
  msg << "(";
  switch (state.type) {
    case ThreadType::Master:
      msg << "master";
      break;
    case ThreadType::ExclusivelyOwned:
      msg << "exclusively-owned";
      break;
    case ThreadType::NotOwned:
      msg << "not-owned";
      break;
    case ThreadType::Shutdown:
      msg << "shutdown";
      break;
  }
  msg << ",no-of-existing-lock-tasks" << state.numberOfExistingLockTasks << ")";
  return msg.str();
}
