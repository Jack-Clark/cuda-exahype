// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _PEANO_DATA_TRAVERSAL_TASK_SET_H_
#define _PEANO_DATA_TRAVERSAL_TASK_SET_H_


#ifdef SharedTBB
#include <tbb/task.h>
#endif


#include <functional>


namespace peano {
  namespace datatraversal {
    class TaskSet;
  }
}


/**
 * Task Set
 *
 * This is my very simple wrapper of task parallelism in Peano. The class
 * provides only one type of constructor. This constructor is passed a set of
 * instances representing one functor each. The constructor executes all of
 * them in parallel (if this is the user's wish) and returns. Consequently,
 * this wrapper is the right choice iff all the task are totally independent
 * of each other.
 *
 * The standard variants do not copy the task you hand in and I expect the
 * functors to be handed in as references. If the tasks are not independent
 * of each other, they thus might induce data races.
 *
 * Please see the documentation of the constructor with only one argument for
 * further details, as this one behaves differently.
 *
 * @author Tobias Weinzierl
 */
class peano::datatraversal::TaskSet {
  private:
    #ifdef SharedTBB
    static tbb::task_group_context  _backgroundTaskContext;


    template <class Functor>
    class GenericTaskWithCopy: public tbb::task {
      private:
        /**
         * See the outer class description for an explanation why this is an
         * attribute, i.e. why we copy the functor here always.
         */
        Functor   _functor;
      public:
        GenericTaskWithCopy(const Functor& functor);
        tbb::task* execute();
    };
    #endif

  public:
    /**
     * Spawn One Asynchronous Task
     *
     * Different to other tasks, I have to copy the functor here. Otherwise,
     * this operation might return, the calling code might destroy the functor,
     * and the asynchronous task then tries to invoke it. This would result in
     * a seg fault.
     *
     * As a consequence, you have to very carefully when the destructor is not
     * empty. Task objects might be copied multiple times and you never know
     * which destructor is the one belonging to the task object that is really
     * executed.
     *
     * <h2> TBB </h2>
     *
     * I do not use the spawn command of TBB here but the enqueue. In
     * particular on Windows systems, I often encounter starvation processes if
     * a load vertices or store vertices task splits very often and spawns too
     * many task children due to this operation.
     *
     * <h2> Deadlocks </h2>
     *
     * Please note that your code might deadlock if you spawn a task without
     * multicore support and if you hope/rely on the fact that this task cannot
     * complete right at the moment but will later on be able to do so.
     */
    template <class Functor>
    inline TaskSet(
      Functor&  task
    );

    /**
     * Invoke operations in parallel. Works fine with lambda
     * expressions:
     *
  peano::datatraversal::TaskSet runParallelTasks(
    [&]() -> void {

    ...
    },
    [&]() -> void {

    ...
    },
    true
  );
     *
     */
    TaskSet(
      std::function<void ()>&& function1,
      std::function<void ()>&& function2,
      bool                     parallelise
    );

    TaskSet(
      std::function<void ()>&& function1,
      std::function<void ()>&& function2,
      std::function<void ()>&& function3,
      bool                     parallelise
    );

    TaskSet(
      std::function<void ()>&& function1,
      std::function<void ()>&& function2,
      std::function<void ()>&& function3,
      std::function<void ()>&& function4,
      bool                     parallelise
    );

    TaskSet(
      std::function<void ()>&& function1,
      std::function<void ()>&& function2,
      std::function<void ()>&& function3,
      std::function<void ()>&& function4,
      std::function<void ()>&& function5,
      bool                     parallelise
    );
};


#include "peano/datatraversal/TaskSet.cpph"

#endif

