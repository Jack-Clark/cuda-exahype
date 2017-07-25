package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.*;


public class SetupBuildEnvironment extends DepthFirstAdapter {
  public Boolean valid = true;

  private java.io.BufferedWriter _writer;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private boolean _requiresFortran;
  private boolean _useOptimisedKernels = false; //at least one solver uses optimised kernels

  private String _likwidInc;
  private String _likwidLib;
  private String _ipcmInc;
  private String _ipcmLib;

  public SetupBuildEnvironment(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
  }

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    int dimensions = Integer.parseInt( node.getDimension().toString().trim() );
    if (dimensions==2) {
      System.out.print("2d experiment ... ok\n");
      try {
        _writer.write("PROJECT_CFLAGS+=-DDim2\n");
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        exc.printStackTrace();
        valid = false;
      }
    }
    else if (dimensions==3) {
      System.out.print("3d experiment ... ok\n");
      try {
        _writer.write("PROJECT_CFLAGS+=-DDim3\n");
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        exc.printStackTrace();
        valid = false;
      }
    }
    else {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
      valid = false;
    }
  }

  @Override
  public void inASharedMemory(ASharedMemory node) {
    try {
      _writer.write("ifeq ($(SHAREDMEM),)\n");
      _writer.write("  SHAREDMEM=TBB\n");
      _writer.write("endif\n");
      System.out.print("shared memory ... TBB (switch to OpenMP manually as indicated below)\n");

      if (!System.getenv().containsKey("TBB_INC")) {
        System.out.print(
            "WARNING: environment variable TBB_INC not set but required if code is built with TBB\n");
      }
      if (!System.getenv().containsKey("TBB_SHLIB")) {
        System.out.print(
            "WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB\n");
      }
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }

  @Override
  public void inADistributedMemory(ADistributedMemory node) {
    try {
      _writer.write("ifeq ($(DISTRIBUTEDMEM),)\n");
      _writer.write("  DISTRIBUTEDMEM=MPI\n");
      _writer.write("endif\n");
      System.out.print("mpi ... switched on \n");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }

  @Override
  public void inAProject(AProject node) {
    _requiresFortran = false;

    try {
      java.io.File logFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/Makefile");

      _writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));

      _writer.write("# *********************************************************************************************"      + "\n");
      _writer.write("# README"                                                                                             + "\n");
      _writer.write("# *********************************************************************************************"      + "\n");
      _writer.write("#"                                                                                                    + "\n");
      _writer.write("# Available configuration Parameters for ExaHyPE"                                                     + "\n");
      _writer.write("#"                                                                                                    + "\n");
      _writer.write("# export variable  |  default-value  |  further values         |  description"                        + "\n");
      _writer.write("#--------------------------------------------------------------------------------------------------"  + "\n");
      _writer.write("# ARCHITECTURE        CPU               Phi, KNL, HSW             Hardware-platform"                 + "\n");
      _writer.write("# COMPILER            Intel             GNU                       Used compiler (and linker)"         + "\n");
      _writer.write("# MODE                Release           Debug, Profile, Asserts   Verbosity and Debug level"          + "\n");
      _writer.write("# SHAREDMEM           None              OMP, TBB                  Shared-memory parallelisation"      + "\n");
      _writer.write("# DISTRIBUTEDMEM      None              MPI                       Distributed-memory parallelisation" + "\n");
      _writer.write("# BOUNDARYCONDITIONS  None              Periodic                  Type of boundary conditions"        + "\n");
      _writer.write("# *********************************************************************************************"      + "\n");

      _writer.write("PEANO_KERNEL_PEANO_PATH=" + _directoryAndPathChecker.peanoKernelPath.getAbsolutePath() + "/peano\n");
      _writer.write("PEANO_KERNEL_TARCH_PATH=" + _directoryAndPathChecker.peanoKernelPath.getAbsolutePath() + "/tarch\n");
      _writer.write("PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=" + _directoryAndPathChecker.peanoToolboxPath.getAbsolutePath() + "/multiscalelinkedcell\n");
      _writer.write("PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=" + _directoryAndPathChecker.peanoToolboxPath.getAbsolutePath() + "/sharedmemoryoracles\n");
      _writer.write("PEANO_TOOLBOX_MPI_BLANCING_PATH=" + _directoryAndPathChecker.peanoToolboxPath.getAbsolutePath() + "/mpibalancing\n");

      _writer.write(
          "EXAHYPE_PATH=" + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "\n");
      _writer.write(
          "PROJECT_PATH=" + _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "\n");
      _writer.write("EXECUTABLE=ExaHyPE-" + node.getName() + "\n");
      _writer.write("\n\n");

      String architecture = "noarch";
      if (node.getArchitecture()!=null) {
        architecture = node.getArchitecture().toString().trim().toLowerCase();
      }

      if (architecture.equals("wsm")) {
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=16");
      } else if (architecture.equals("snb")) {
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=32");
      } else if (architecture.equals("hsw")) {
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=32");
      } else if (architecture.equals("knc")) {
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=64");
      } else if (architecture.equals("knl")) {
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=64");
      } else {
        // noarch or unsupported architecture or undefined
        _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=16");
      }

      _writer.write("\n");
      _writer.write("# Several MPI variants face problems with multithreaded MPI. As we run into \n");
      _writer.write("# such issues multiple times, we disable by default multithreaded MPI in ExaHyE. \n");
      _writer.write("# However, feel free to give it a try in your code on your system by disabling \n");
      _writer.write("# this flag. \n");
      _writer.write("PROJECT_CFLAGS+=-DnoMultipleThreadsMayTriggerMPICalls\n");

    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    if (node.getLanguage().getText().trim().equals("C")) {
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      _requiresFortran = true;
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
    }
    
    _useOptimisedKernels = _useOptimisedKernels 
                            || (node.getLanguage().getText().trim().equals("C") 
                                && (node.getKernel().toString().trim().equals(eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier)
                                    ||  node.getKernel().toString().trim().equals(eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier)));
    
  }

  @Override
  public void inAProfiling(AProfiling node) {
    if (node.getLikwidInc() != null) {
      _likwidInc = node.getLikwidInc().toString().trim();
    }

    if(node.getLikwidLib() != null) {
      _likwidLib = node.getLikwidLib().toString().trim();
    }

    if (node.getIpcmInc() != null) {
      _ipcmInc = node.getIpcmInc().toString().trim();
    }
    if (node.getIpcmLib() != null) {
      _ipcmLib = node.getIpcmLib().toString().trim();
    }
  };

  @Override
  public void outAProject(AProject node) {

    try {
      _writer.write("\n\n");
      if (_requiresFortran) {
        _writer.write("ifeq ($(MIXEDLANG),)\n");
        _writer.write("  MIXEDLANG=Yes\n");
        _writer.write("endif\n");
      }
      if (_useOptimisedKernels) {
        _writer.write("ifneq ($(call tolower,$(MODE)),release)\n");
        _writer.write(" PROJECT_CFLAGS += -DTEST_OPT_KERNEL\n");
        _writer.write("endif\n");
        _writer.write("\n\n");
        _writer.write("ifeq ($(COMPILE_OPT_KERNEL),)\n");
        _writer.write("  COMPILE_OPT_KERNEL=Yes\n");
        _writer.write("endif\n");
        String architecture = node.getArchitecture().toString().trim().toLowerCase();
        if(!architecture.equals("noarch")) _writer.write("ARCHITECTURE="+architecture+"\n");
      }
      if (_likwidInc != null) {
        _writer.write("PROJECT_CFLAGS+=-I" + _likwidInc + "\n");
        _writer.write("PROJECT_CFLAGS+=-DLIKWID_AVAILABLE\n");
      }
      if (_ipcmInc != null) {
        _writer.write("PROJECT_CFLAGS+=-I" + _ipcmInc + "\n");
        _writer.write("PROJECT_CFLAGS+=-DIPCM_AVAILABLE\n");
      }
      if (_likwidLib != null) {
        _writer.write("PROJECT_LFLAGS+=-L" + _likwidLib + " -llikwid\n");
      }
      if (_ipcmLib != null) {
        _writer.write("PROJECT_LFLAGS+=-L" + _ipcmLib + " -lintelpcm\n");
      }
      _writer.write("\n\n");
      _writer.write(
          "-include " + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "/Makefile\n");
      _writer.write("\n\n\n\n");
      _writer.write("all: \n");
      _writer.write("\t@echo " + node.getName() + "\n");
      _writer.write("\t@echo =================\n");
      _writer.write("\t@echo An ExaHyPE solver\n");

      System.out.print("store pathes and default settings in makefile ... ok");
      System.out.println("\n\n\n\n");
      System.out.print("please change into directory "
          + _directoryAndPathChecker.outputDirectory.getAbsolutePath() + " and type make \n");
      System.out.print("ensure that you set all environment variables before:\n");
      System.out.print("  export COMPILER=GNU  \t\t\tSelect GNU compiler\n");
      System.out.print("  export COMPILER=Intel\t\t\tSelect Intel compiler (default)\n");
      System.out.print("\n");
      System.out.print("  export MODE=Debug\t\t\tBuild debug version of code\n");
      System.out.print(
          "  export MODE=Asserts\t\t\tBuild release version of code that is augmented with assertions\n");
      System.out.print(
          "  export MODE=Profile\t\t\tBuild release version of code that produces profiling information\n");
      System.out.print("  export MODE=Release\t\t\tBuild release version of code (default)\n");
      System.out.print("\n");
      System.out.print(
          "  export SHAREDMEM=TBB\t\t\tUse Intel's Threading Building Blocks (TBB) for shared memory parallelisation\n");
      System.out.print(
          "  export SHAREDMEM=OMP\t\t\tUse OpenMP for shared memory parallelisation\n");
      System.out.print(
          "  export SHAREDMEM=None\t\t\tDo not use shared memory (default if not indicated otherwise by \"shared memory ...\" message above)\n");
      System.out.print("\n");
      System.out.print("  export DISTRIBUTEDMEM=MPI\t\tUse MPI\n");
      System.out.print("  export DISTRIBUTEDMEM=None\t\tDo not use MPI (default)\n");
      System.out.print("\n");
      System.out.print(
          "  export TBB_INC=-I...\t\t\tIndicate where to find TBB headers (only required if SHAREDMEM=TBB). Please add -I (Linux) prefix to path\n");
      System.out.print(
          "  export TBB_SHLIB=\"-L... -ltbb\"\tIndicate where to find TBB's shared libraries (only required if SHAREDMEM=TBB). Variable has to comprise both search path and library name\n");
      System.out.print("\n\n");
      System.out.print(
          "  If SHAREDMEM or DISTRIBUTEDMEM are not specified, they fall back to \"None\".\n");
      System.out.print("\n");
      System.out.print(
          "  If you run CSH, please replace \"export ARG=VALUE\" with \"setenv ARG VALUE\".\n");

      _writer.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      exc.printStackTrace();
      valid = false;
    }
  }
}
