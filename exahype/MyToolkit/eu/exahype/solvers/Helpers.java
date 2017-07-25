package eu.exahype.solvers;

import java.io.IOException;

public class Helpers {
  
  public static void writeMinimalADERDGSolverHeader(
      String solverName, java.io.BufferedWriter writer, String projectName, boolean hasConstants,
      int order, int dimensions, int numberOfUnknowns, int numberOfParameters, boolean enableProfiler) throws IOException {
    writeHeaderCopyright(writer);
    writeHeaderIncludesAndDefines(writer, solverName, projectName);
    writeHeaderMinimalADERDGClassSignature(writer, solverName, projectName, hasConstants, order, dimensions, numberOfUnknowns, numberOfParameters, enableProfiler);
  }

  public static void writeMinimalFiniteVolumesSolverHeader(
	      String solverName, java.io.BufferedWriter writer, String projectName, boolean hasConstants) throws IOException {
    writeHeaderCopyright(writer);
    writeHeaderIncludesAndDefines(writer, solverName, projectName);
    writeHeaderMinimalFiniteVolumesClassSignature(writer, solverName, projectName, hasConstants);
  }

  
  private static void writeHeaderMinimalADERDGClassSignature(
      java.io.BufferedWriter writer, String solverName, String projectName, boolean hasConstants,
      int order, int dimensions, int numberOfUnknowns, int numberOfParameters, boolean enableProfiler) throws IOException {
    writer.write(
        "class " + projectName + "::" + solverName + ": public exahype::solvers::ADERDGSolver {\n");
    writer.write("  public:\n");
    writer.write("\n");
    writer.write("    // Sorry for being inconsistent here: While AderDGSolver offers the methods getNumberOfVariables() etc.,\n");
    writer.write("    // in static context they cannot be accessed. Thus the toolkit offers you access to the variables here.\n");
    writer.write("    // Thank you, Toolkit!\n");
    writer.write("    static const int nVar = "+numberOfUnknowns+";\n");
    writer.write("    static const int nDim = "+dimensions+";\n");
    writer.write("    static const int nParams = "+numberOfParameters+";\n");
    writer.write("    static const int order = "+order+";\n");
    writer.write("\n");
    
    
    if (hasConstants) {
      writer.write("    " + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping"+
    (enableProfiler ? ", std::unique_ptr<exahype::profilers::Profiler> profiler" : "")+
    ", std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView constants);\n\n");
    }
    else {
      writer.write("    " + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping"+
    (enableProfiler ? ", std::unique_ptr<exahype::profilers::Profiler> profiler" : "")+
    ", std::vector<std::string>& cmdlineargs);\n\n");
    }

    writer.write(
        "    void spaceTimePredictor(double* lQhbnd,double* lFhbnd, double** tempSpaceTimeUnknowns, double** tempSpaceTimeFluxUnknowns, double*  tempUnknowns, double*  tempFluxUnknowns, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt) override; \n");
    writer.write(
        "    void solutionUpdate(double* luh, const double* const lduh, const double dt) override;\n");
    writer.write(
        "    void volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) override;\n");
    writer.write(
        "    void surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) override;\n");
    writer.write(
        "    void riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) override;\n");
    writer.write(
        "    void boundaryConditions(double* fluxOut, double* stateOut, const double* const fluxIn, const double* const stateIn, const tarch::la::Vector<DIMENSIONS, double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>& cellSize, const double t,const double dt, const int faceIndex, const int normalNonZero) override;\n");
    writer.write(
        "    double stableTimeStepSize(const double* const luh, double* tempEigenvalues, const tarch::la::Vector<DIMENSIONS,double>& dx) override;\n");
    writer.write(
        "    void solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;\n");
    writer.write(
        "    bool hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) override;\n");
    writer.write(
        "    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) override;\n");
    writer.write(
        "    void faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;\n");
    writer.write(
        "    void faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) override;\n");
    writer.write(
        "    void volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;\n");
    writer.write(
        "    void volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) override;\n");
  }

  /**
   * Creates all the public operations that are mandatory for any solver.
   */
  private static void writeHeaderMinimalFiniteVolumesClassSignature(
      java.io.BufferedWriter writer, String solverName, String projectName, boolean hasConstants) throws IOException {
    writer.write(
        "class " + projectName + "::" + solverName + ": public exahype::solvers::FiniteVolumesSolver {\n");
    writer.write("  public:\n");
    writer.write("    " + solverName + "(int cellsPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::vector<std::string>& cmdlineargs, std::unique_ptr<exahype::profilers::Profiler> profiler);\n\n");

    writer.write("    double stableTimeStepSize( double* luh[THREE_POWER_D], const tarch::la::Vector<DIMENSIONS, double>& dx) override; \n\n" );
    writer.write("    void   solutionAdjustment( double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) override; \n\n");
    writer.write("    bool   hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t) override; \n\n" );
    writer.write("    exahype::solvers::Solver::RefinementControl refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) override; \n\n" );
    writer.write("    void solutionUpdate(double* luh[THREE_POWER_D], const tarch::la::Vector<DIMENSIONS, double>& dx, const double dt, double& maxAdmissibleDt) override; \n\n" );
  }

  /**
   * Write header with ExaHyPE copyright. Should be inserted for any solver's
   * header.
   */
  public static void writeHeaderCopyright(java.io.BufferedWriter writer) throws IOException {
    writer.write("// This file is generated by the ExaHyPE toolkit.\n");
    writer.write("// Please do not modify - it will be overwritten by the next\n");
    writer.write("// ExaHyPE toolkit call.\n");
    writer.write("// \n");
    writer.write("// ========================\n");
    writer.write("//   www.exahype.eu\n");
    writer.write("// ========================\n");
  }

  /**
   * Adds all the default includes of any solver as well as the solver define.
   * Is used by all solvers.
   */
  private static void writeHeaderIncludesAndDefines(
      java.io.BufferedWriter writer, String solverName, String projectName) throws IOException {
    writer.write("\n\n");
    writer.write("#include <memory>\n\n");
    writer.write("#include \"exahype/Parser.h\"\n");
    writer.write("#include \"exahype/profilers/Profiler.h\"\n");
    writer.write("#include \"exahype/solvers/ADERDGSolver.h\"\n");
    writer.write("#include \"exahype/solvers/FiniteVolumesSolver.h\"\n");
    writer.write("\n\n\n");

    writer.write("namespace " + projectName + "{\n");
    writer.write("  class " + solverName + ";\n");
    writer.write("}\n\n\n");
  }

  public static void writeMinimalADERDGSolverUserImplementation(String solverName,
      java.io.BufferedWriter writer, String projectName, int numberOfVariables, int numberOfParameters, int order, boolean hasConstants)
      throws IOException {
    writer.write("#include \"" + solverName + ".h\"\n\n");
    writer.write("#include <memory>\n\n");

    // init
    writer.write("void " + projectName + "::" + solverName +
      "::init(std::vector<std::string>& cmdlineargs"+(hasConstants?", exahype::Parser::ParserView& constants":"")+") {\n");
    writer.write("  // This function is called by the constructor.\n");
    writer.write("  // You can access spec file parameters as well as command line arguments (argv as std::vector).\n");
    writer.write("  // @todo Please implement/augment if required.\n");
    writer.write("}\n\n");
  }
  
  public static void writeMinimalFiniteVolumesSolverUserImplementation(String solverName,
      java.io.BufferedWriter writer, String projectName, int numberOfVariables, int numberOfParameters, int patchSize, boolean hasConstants)
      throws IOException {
    writer.write("#include \"" + solverName + ".h\"\n\n");
    writer.write("#include <memory>\n\n");

    // init
    writer.write("void " + projectName + "::" + solverName +
      "::init() {\n");
    writer.write("  // This function is called inside the constructur.\n");
    writer.write("  // @todo Please implement/augment if required.\n");
    writer.write("}\n\n");
  }

  static public void invokeCodeGenerator(String solverName, int numberOfUnknowns, int numberOfParameters, int order,
      boolean isLinear, int dimensions, String microarchitecture, String pathToLibxsmm)
      throws IOException {
    String currentDirectory = System.getProperty("user.dir");
    java.io.File pathToCodeGenerator =
        new java.io.File(currentDirectory + "/CodeGenerator/Driver.py");
    if (!pathToCodeGenerator.exists()) {
      System.err.println("ERROR: Code generator not found. Can't generate optimised kernels. Path: " + pathToCodeGenerator.toString());
      throw new IOException();
    }
    
    if(pathToLibxsmm == null || pathToLibxsmm.isEmpty()) {
      System.err.println("ERROR: Libxsmm path not specified");
      throw new IOException();
    }
    
    java.io.File pathToLibxsmmMakefile = //To test if the libxsmm folder is correct
        new java.io.File(java.nio.file.Paths.get(currentDirectory,pathToLibxsmm,"Makefile").toString());
    if (!pathToLibxsmmMakefile.exists()) {
      System.err.println("ERROR: Libxsmm makefile not found. Can't generate optimised kernels. Path: " + pathToLibxsmmMakefile.toString());
      throw new IOException();
    }

    String numericsParameter = isLinear ? "linear" : "nonlinear";

    // set up the command to execute the code generator
    String args = " " + "Euler" + " " + numberOfUnknowns + " " + order + " " //TODO JMG see why Euler instead of solverName
        + Integer.toString(dimensions) + " " + numericsParameter + " " + microarchitecture + " "
        + currentDirectory + "/"  + pathToLibxsmm + " "; 

    String bashCommand = "env python3 " + pathToCodeGenerator + args;

    Runtime runtime = Runtime.getRuntime();
	System.out.println("Codegenerator command line: "+bashCommand);
    // execute the command line program
    Process codeGenerator = runtime.exec(bashCommand);

    // capture any output that is produced by the code generator and print it line-by-line
    java.io.InputStream stdout = codeGenerator.getInputStream();
    java.io.BufferedReader stdoutReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stdout));
    String line = "";
    while ((line = stdoutReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }
    java.io.InputStream stderr = codeGenerator.getErrorStream();
    java.io.BufferedReader stderrReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stderr));
    while ((line = stderrReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }

    // in order to stop further toolkit execution if the code generator fails,
    // explicitly wait for the process
    try {
        int exitValue = codeGenerator.waitFor();
        if(exitValue != 0) {
            System.err.println("ERROR: Code Generator failed with exit value " + exitValue);
            throw new IOException(); // <- also done in line 186. This is abusing the exception system.
        }
    } catch(InterruptedException e) {
        System.err.println("This is very bad. I don't know what's going on.");
        throw new IOException();
    }
  } // invokeCodeGenerator
}
