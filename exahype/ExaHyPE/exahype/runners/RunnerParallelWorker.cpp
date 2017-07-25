/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#include "exahype/runners/Runner.h"

#ifdef Parallel
#include "exahype/repositories/Repository.h"
#include "peano/parallel/messages/ForkMessage.h"
#include "peano/utils/Globals.h"
#include "peano/utils/UserInterface.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/parallel/NodePool.h"

int exahype::runners::Runner::runAsWorker(
    exahype::repositories::Repository& repository) {
  int newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  while (newMasterNode !=
         tarch::parallel::NodePool::JobRequestMessageAnswerValues::Terminate) {
    if (newMasterNode >=
        tarch::parallel::NodePool::JobRequestMessageAnswerValues::NewMaster) {
      peano::parallel::messages::ForkMessage forkMessage;
      forkMessage.receive(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          tarch::parallel::NodePool::getInstance().getTagForForkMessages(),
          true, ReceiveIterationControlMessagesBlocking);

      repository.restart(
          forkMessage.getH(), forkMessage.getDomainOffset(),
          forkMessage.getLevel(),
          forkMessage.getPositionOfFineGridCellRelativeToCoarseGridCell());

      bool continueToIterate = true;
      while (continueToIterate) {
        switch (repository.continueToIterate()) {
          case exahype::repositories::Repository::Continue:
            {
              repository.iterate();
              logInfo("startNewTimeStep(...)",
                "\tmemoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
            }
            break;
          case exahype::repositories::Repository::Terminate:
            continueToIterate = false;
            break;
          case exahype::repositories::Repository::RunGlobalStep:
            runGlobalStep();
            break;
        }
      }

      // insert your postprocessing here
      // -------------------------------

      // -------------------------------

      repository.terminate();
    } else if (newMasterNode ==
               tarch::parallel::NodePool::JobRequestMessageAnswerValues::
                   RunAllNodes) {
      runGlobalStep();
    }
    newMasterNode = tarch::parallel::NodePool::getInstance().waitForJob();
  }
  return 0;
}

void exahype::runners::Runner::runGlobalStep() {
  // You might want to remove this assertion, but please consult the
  // documentation before.
  assertion(!peano::parallel::loadbalancing::Oracle::getInstance()
                 .isLoadBalancingActivated());

  // insert yourcode here
  // -------------------------------

  // -------------------------------
}
#endif
