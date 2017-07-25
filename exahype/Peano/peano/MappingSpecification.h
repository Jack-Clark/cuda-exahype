// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _PEANO_MAPPING_SPECIFICATION_H_
#define _PEANO_MAPPING_SPECIFICATION_H_


#include <string>


namespace peano {
  struct MappingSpecification;
}


/**
 * Specification for a mapping operation
 *
 * Gives Peano hints on the behaviour and, hence, supports a more aggressive
 * on-the-fly optimisation.
 *
 *
|| manipulates || Nop | Operation is empty, and the algorithm's semantics is preserved, if Peano kicks out the whole function call.
||             || OnlyLeaves | Algorithm's semantics is preserved, if Peano kicks out the whole function call for refined vertices of cells that have solely refined vertices around them.
||             || WholeTree  | No events are eliminated.
|| multithreading || Serial                    | Do not run these operation in parallel on a shared memory machine.
||                ||                           | Application hence doesn't need semaphores.
||                || RunConcurrentlyOnFineGrid | Peano runs the events in parallel on the finest grid, i.e. doesn't
||                ||                           | care about any data dependencies. For multiscale events such as
||                ||                           | ascend, this operation runs the code in parallel for all coarse
||                ||                           | grid cells.
||                || AvoidCoarseGridRaces      | Peano can try to speed up the application due to multithreading.
||                ||                           | It ensures that events are invoked such that the coarse grid data
||                ||                           | (of enterCell, e.g.) is not shared with another thread. Is more
||                ||                           | restrictive than AvoidFineGridRaces, i.e. ensures this data
||                ||                           | consistency as well. Usually leads to colouring with @f$ 6^d @f$ or
||                ||                           | @f$ 7^d @f$ colors.
||                || AvoidFineGridRaces        | Peano can try to speed up the application due to multithreading.
||                ||                           | However, Peano ensures that events are invoked such that the fine
||                ||                           | grid data (of enterCell, e.g.) is not shared with another thread.
||                ||                           | Introduces @f$ 2^d @f$ colouring on the cells.
||                ||                           | Variant is undefined (and thus may not be chosen) for inter-level
||                ||                           | events such as ascend.
||                ||                           | Variant seems to be unnecessary for touchVertex...Time events as
||                ||                           | those are not passed their neighbouring vertices anyway, i.e. there
||                ||                           | may not be any data races. However, if you work with pointers to
||                ||                           | adjacent cells, it can make sense to have this colouring.
 *
 * The order is Serial>AvoidCoarseGridRaces>AvoidFineGridRaces>RunConcurrentlyOnFineGrid.
 * If two mappings are combined one holding AvoidCoarseGridRaces and one holding
 * RunConcurrentlyOnFineGrid, the combination holds AvoidCoarseGridRaces.
 *
 * The state boolean is actually only required if you run the code with
 * shared memory parallelisation. On a multicore chip, the kernel may invoke
 * events in parallel. If a mapping's operation modifies the state, the code
 * has to copy the mapping to each thread before it starts the parallel
 * invocation, and it has to reduce the mappings via mergeWithWorker(...)
 * in the end. If the operations do not alter the state, we can skip the
 * copying and the reduction. We basically may assume that the event is
 * const.
 *
 * @author Tobias Weinzierl
 */
struct peano::MappingSpecification {
  enum Manipulates  {
    Nop, OnlyLeaves, WholeTree
  };

  enum Multithreading {
    Serial, AvoidCoarseGridRaces, AvoidFineGridRaces, RunConcurrentlyOnFineGrid
  };

  Manipulates     manipulates;
  Multithreading  multithreading;
  bool            altersState;

  MappingSpecification(Manipulates manipulates_, Multithreading multithreading_, bool altersState_);

  /**
   * Most general specification
   *
   * Is used by the adapters to merge multiple specifications.
   */
  static MappingSpecification getMinimalSpecification();

  std::string toString() const;
};

/**
 * Combine two specifications. The resulting specification has weaker or equal
 * constraints. If one specification work for example only on the leaves and
 * supports multicore parallelism and the other one works on the whole tree
 * serially, then the result specification works serially on the whole tree.
 */
peano::MappingSpecification operator&(const peano::MappingSpecification& lhs, const peano::MappingSpecification& rhs);


#endif
