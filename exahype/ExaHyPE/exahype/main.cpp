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
 
#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"

#include "exahype/Parser.h"
#include "exahype/runners/Runner.h"
#include "buildinfo.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include <vector>
#include <string>
#include <cstdlib> // getenv
#include <iostream>
#include <cstdio>

tarch::logging::Log _log("");

void version(); // version dumping, see below
void help(const char* programname);  // A help message

int main(int argc, char** argv) {
  peano::fillLookupTables();

  //
  //   Setup environment
  // =====================
  //
  int parallelSetup = peano::initParallelEnvironment(&argc, &argv);
  if (parallelSetup != 0) {
#ifdef Parallel
    // Please do not use the logging if MPI doesn't work properly.
    std::cerr << "mpi initialisation wasn't successful. Application shut down"
              << std::endl;
#else
    _log.error("main()",
               "mpi initialisation wasn't successful. Application shut down");
#endif
    return parallelSetup;
  }

  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
  if (sharedMemorySetup != 0) {
    logError("main()",
             "shared memory initialisation wasn't successful. Application shut "
             "down");
    return sharedMemorySetup;
  }

  //
  //   Parse config file
  // =====================
  //
  if (argc < 2) {
    logError("main()", "Usage: ./ExaHyPE config-file [additional args passed to Solver...]");
    return -1;
  }

  if(std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
    help(argv[0]);
    return -1;
  }

  // TODO(Dominic): This is not a reboust way to parse in the command line arguments; use the vector below instead.
  bool onlyShowVersion=std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v";

  exahype::Parser parser;
  parser.readFile(argv[onlyShowVersion ? 2 : 1]);

  if (!parser.isValid()) {
    logError("main()", "invalid config file. Quit");
    return -2;
  }
  
  // Collect all command line arguments for the Solvers
  std::vector<std::string> cmdlineargs(argv + 1, argv + argc);

  //
  //   Init solver registries
  // =====================================
  //
  kernels::initSolvers(parser, cmdlineargs);

  if (onlyShowVersion) {
    version();
    kernels::finalise();
    return 0;
  }

  //
  //   Configure the logging
  // =========================
  //
  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  #ifdef Parallel
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      true,   // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #elif defined(Asserts) || defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      true,   // logTrace
      parser.getLogFileName() );
  #else
  tarch::logging::CommandLineLogger::getInstance().setLogFormat(
      " ",    // columnSeparator
      true,   // logTimeStamp
      false,  // logTimeStampHumanReadable
      false,  // logMachineName
      true,   // logMessageType
      false,   // logTrace
      parser.getLogFileName() );
  #endif

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  if (!tarch::logging::LogFilterFileReader::parsePlainTextFile(
          "exahype.log-filter")) {
    tarch::logging::CommandLineLogger::getInstance().clearFilterList();
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("info", false));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry(
            "info", -1, "peano::grid", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", true));
    tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
        ::tarch::logging::CommandLineLogger::FilterListEntry("debug", -1,
                                                             "exahype", false));
  }

//
//   Run tests
// =============
//
#if defined(Debug) || defined(Asserts)
if(! std::getenv("EXAHYPE_SKIP_TESTS")) { // cf issue #74
  tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
  int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
                         .getTestCaseCollection()
                         .getNumberOfErrors();

  if (testExitCode != 0) {
    logError("main()", "unit tests failed. Quit.");
    return -2;
  }
} else {
  logInfo("main()", "Skipping tests as EXAHYPE_SKIP_TESTS is set."
     "We do so because tests are broken in the moment and nobody repairs them."); //  TODO(Sven,Dominic,JM): Fix tests.
} // end if getenv(EXAHYPE_SKIP_TESTS)
#endif

  exahype::runners::Runner runner(parser);
  int programExitCode = runner.run();

  if (programExitCode == 0) {
#ifdef Parallel
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logInfo("main()", "Peano terminates successfully");
    }
#else
    logInfo("main()", "Peano terminates successfully");
#endif
  } else {
    logInfo("main()", "quit with error code " << programExitCode);
  }

  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();
  peano::releaseCachedData();

  kernels::finalise();

  return programExitCode;
}


void version() {
  std::cout << "This is an ExaHyPE executable (http://exahype.eu)\n";
  std::cout << "Compiled on host " << EXAHYPE_BUILD_HOST << " at " << EXAHYPE_BUILD_DATE << "\n";
#ifdef EXAHYPE_GIT_INFO
  std::cout << "ExaHyPE git version: " << EXAHYPE_GIT_INFO << "\n";
#else
  std::cout << "ExaHyPE git version: n/a\n";
#endif
#ifdef PEANO_SVN_INFO
  std::cout << "Peano svn version:   " << PEANO_SVN_INFO << "\n";
#else
  std::cout << "Peano svn version:   n/a\n";
#endif
  std::cout << "\n";

  std::cout << "Compile time options\n";
  std::cout << "====================\n";
#ifdef DIMENSIONS
  std::cout << "Dimensions:    "<< DIMENSIONS << "\n";
#else
  std::cout << "Dimensions:    not determinable!\n";
#endif

#ifdef Debug
  std::cout << "Debug:         YES\n";
#else
  std::cout << "Debug:         no\n";
#endif
  
#ifdef Asserts
  std::cout << "Assertions:    YES\n";
#else
  std::cout << "Assertions:    no\n";
#endif

#ifdef Parallel
  std::cout << "MPI Support:   YES\n";
#else
  std::cout << "MPI Support:   no\n";
#endif
  
#ifdef EXAHYPE_CFL_FACTOR // issue #100
  std::cout << "CFL Factor:    "<< EXAHYPE_CFL_FACTOR << "\n";
#else
  std::cout << "CFL Factor:    Default (0.9 or so)\n";
#endif

  std::cout << "\n";
  std::cout << "Makesystem build options\n";
  std::cout << "========================\n";
#ifdef EXAHYPE_BUILDINFO_AVAILABLE
  std::cout << EXAHYPE_BUILD_INFO << "\n";
#else
  std::cout << "Symbols n/a" << "\n";
#endif

  std::cout << "\n";
  std::cout << "Toolkit static registry info\n";
  std::cout << "============================\n";
  kernels::toString(std::cout);
  std::cout << "\n";
}

void help(const char* programname) {
  std::cout << "Usage: " << programname << " <YourApplication.exahype>\n";
  std::cout << "\n";
  std::cout << "   where YourApplication.exahype is an ExaHyPE specification file.\n";
  std::cout << "   Note that you should have compiled ExaHyPE with this file as there\n";
  std::cout << "   are some compile time constants.\n";
  std::cout << "\n";
  std::cout << "   Other possible parameters:\n";
  std::cout << "\n";
  std::cout << "    --help     Show this help message\n";
  std::cout << "    --version  Show version and other hard coded information\n";
  std::cout << "\n";
}
