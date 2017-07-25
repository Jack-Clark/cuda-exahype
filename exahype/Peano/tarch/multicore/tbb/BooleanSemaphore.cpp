#if defined(SharedTBB) || defined(SharedTBBInvade)
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/logging/Log.h"


#include <limits>


#include <tbb/tbb_machine.h>
#include <tbb/task.h>
#include <tbb/tbb_thread.h>


tarch::multicore::BooleanSemaphore::BooleanSemaphore() {
}


tarch::multicore::BooleanSemaphore::~BooleanSemaphore() {
}


void tarch::multicore::BooleanSemaphore::enterCriticalSection() {
  _mutex.lock();
}


void tarch::multicore::BooleanSemaphore::leaveCriticalSection() {
  _mutex.unlock();
}


void tarch::multicore::BooleanSemaphore::sendTaskToBack() {
  tbb::this_tbb_thread::yield();
}
#endif
