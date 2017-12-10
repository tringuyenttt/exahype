// @todo Martin

#ifndef _SHMINVADE_SHMCONTROLLER_H_
#define _SHMINVADE_SHMCONTROLLER_H_


#include <tbb/atomic.h>
#include <tbb/spin_mutex.h>
#include <tbb/concurrent_hash_map.h>


// @todo Remove Namespace
namespace shminvade {
  class SHMController;
  class SHMLockTask;
  class SHMStrategy;
  class SHMOccupyAllCoresStrategy;
  class SHMMultipleRanksPerNodeStrategy;
}



class shminvade::SHMController {
  public:
    enum class ThreadType {
      Master,
      ExclusivelyOwned,
      NotOwned,
      Shutdown
    };

    struct ThreadState {
      ThreadState(ThreadType type_):
        type(type_),
        numberOfExistingLockTasks(0) {}

      ThreadType type;

      /**
       * It is really important to ensure that not too many lock threads are
       * flying around in the system. So we count them per thread and do issue
       * lock threads if and only if there is a need to do so. Could obviously
       * be a map over bools, but I just prefer to count up and down today.
       */
      int        numberOfExistingLockTasks;
    };

    static std::string toString( ThreadState state );

  private:
    tbb::atomic<bool> _switchedOn;


/*
    https://software.intel.com/en-us/node/506171
    Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.

    @todo Wir koennen da nicht drueber iterieren, weil sich die States ja staendig aendern koennen
*/
    typedef tbb::concurrent_hash_map<pid_t, ThreadState> ThreadTable;
    ThreadTable  _threads;

    void setThreadTableEntry( pid_t pid, ThreadState state );
    ThreadState getThreadTableEntry( pid_t pid ) const;

    /**
     * Just register the master thread through registerNewThread(). This
     * implies that a lock task is launched for the master thread which
     * eventually will go down - the latest when we shut down the system.
     * Before it does so, it will identify lots of other threads that would
     * in theory be available.
     */
    SHMController();

    /**
     * Insert a new thread entry and then launch the lock task for it
     * immediately.
     *
     * May not be called if thread already is known
     */
    void registerNewThread(int threadId, ThreadType initialstate=ThreadType::NotOwned );

    /**
     * Is typically used by the strategies when they go down.
     */
    void retreatFromAllCores();

    friend class SHMLockTask;
    friend class SHMStrategy;
    friend class SHMOccupyAllCoresStrategy;
    friend class SHMMultipleRanksPerNodeStrategy;
  public:
    ~SHMController();

    static SHMController&  getInstance();

    void switchOn();
    void switchOff();

    int getMaxAvailableCores() const;
    int getFreeCores() const;
    int getBookedCores() const;

    /**
     * The shutdown sets all internal thread states to Shutdown such that all
     * lock threads hanging around in the system shut down, too. It then goes
     * to sleep to allow the lock threads to see this new state being a signal
     * to them. Before, we switch the invasion overall off and thus can then
     * terminate.
     */
    void shutdown();

    void retreatFromCore( pid_t pid );
};

#endif
