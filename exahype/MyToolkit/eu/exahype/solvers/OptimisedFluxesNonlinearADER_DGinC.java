package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class OptimisedFluxesNonlinearADER_DGinC implements Solver {
  public static final String Identifier = "optimised::fluxes::nonlinear";

  private int _dimensions;
  private int _numberOfUnknowns;
  private int _numberOfParameters;
  private int _order;
  private String _microarchitecture;
  private String _pathToLibxsmm;
  private boolean _enableProfiler;
  private boolean _hasConstants;

  public OptimisedFluxesNonlinearADER_DGinC(int dimensions, int numberOfUnknowns, int numberOfParameters, int order,
      String microarchitecture, String pathToLibxsmm, boolean enableProfiler, boolean hasConstants) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
    _microarchitecture = microarchitecture;
    _pathToLibxsmm = pathToLibxsmm;
    _enableProfiler = enableProfiler;
    _hasConstants   = hasConstants;
  }

  @Override
  public void writeAbstractHeader(BufferedWriter writer, String solverName,
      String projectName) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer,
      String solverName, String projectName) throws IOException {
    // TODO Auto-generated method stub
    
  }
  
  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    IncludeOnceHelper ifndef = new IncludeOnceHelper(writer, solverName+"_CLASS_HEADER");
    ifndef.open();
    Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName, _hasConstants, _order, _dimensions, _numberOfUnknowns, _numberOfParameters, 
        _enableProfiler);
    writer.write("  private:\n");
	writer.write(
        "    static void flux(const double* const Q, double** F);\n"); //TODO JMG Remove fluxSplitted when not needed anymore
    if (_dimensions == 2) {
      writer.write("    static void fluxSplitted(const double* const Q, double* f, double* g);\n");
    } else {
      writer.write(
          "    static void fluxSplitted(const double* const Q, double* f, double* g, double* h);\n");
    }
	
    writer.write(
        "    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write(
        "    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");
		
    writer.write("    void init(std::vector<std::string>& cmdlineargs"+(_hasConstants ? ", exahype::Parser::ParserView& constants" : "")+");\n");
    writer.write("    void source(const double* const Q, double* S);\n");
    writer.write("    void boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut);\n");
    writer.write(
        "    void ncp(const double* const Q, const double* const gradQ, double* BgradQ);\n");
    writer.write(
        "    void matrixb(const double* const Q, const int normalNonZero, double* Bn);\n");

    writer.write("};\n\n\n");
	ifndef.close();
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.invokeCodeGenerator(solverName, _numberOfUnknowns, _numberOfParameters, _order, false, _dimensions,
        _microarchitecture, _pathToLibxsmm);

    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/optimised/Kernels.h\"\n");
    writer.write("\n\n\n");
	
	// constructor
    if (_hasConstants) {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView constants):\n");
    }
    else {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::vector<std::string>& cmdlineargs):\n");
    }


    writer.write("  exahype::solvers::ADERDGSolver("
        + "\""+solverName+"\", nVar /* numberOfUnknowns */, "
        + "nParams /* numberOfParameters */, order + 1 "
        + " /* nodesPerCoordinateAxis */, maximumMeshSize, timeStepping) {\n");
    if(_hasConstants) {
       writer.write("  init(cmdlineargs, constants);\n");
    } else {
       writer.write("  init(cmdlineargs);\n");
       writer.write("  // PS: If you miss access to user constants here, enable them in the toolkit\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
	
	//TODO JMG Remove fluxSplitted when not needed anymore
    writer.write("void " + projectName + "::" + solverName
        + "::fluxSplitted(const double* const Q, double* f, double* g"+(_dimensions == 2 ? "" : ", double* h")+") {\n");
    writer.write("   double* F["+_dimensions+"];\n");
	writer.write("   F[0] = f;\n");
	writer.write("   F[1] = g;\n");
	if(_dimensions == 3) 
		writer.write("   F[2] = h;\n");
	writer.write("   flux(Q,F);\n");
    writer.write("}\n");
	writer.write("\n\n\n");
	
    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"spaceTimePredictor\");\n");
    }
    writer.write("   kernels::aderdg::optimised::picardLoop<fluxSplitted>( tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0], luh, dx, dt );\n"); //TODO remove fluxSplitted for flux
    writer.write("   kernels::aderdg::optimised::predictor( tempUnknowns, tempFluxUnknowns, tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0] );\n");
    writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, tempUnknowns, tempFluxUnknowns );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"spaceTimePredictor\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"solutionUpdate\");\n");
    }
    writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"solutionUpdate\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"volumeIntegral\");\n");
    }
    writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"volumeIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"surfaceIntegral\");\n");
    }
    writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"surfaceIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"riemannSolver\");\n");
    }
    writer.write(
        "   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"riemannSolver\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize( const double* const luh, double* tempEigenvalues, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"stableTimeStepSize\");\n");
    }
    writer.write(
        "   double d = kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"stableTimeStepSize\");\n");
    }
    writer.write("   return d;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"solutionAdjustment\");\n");
    }
    writer.write(
        "   kernels::aderdg::optimised::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt );\n");
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"solutionAdjustment\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"faceUnknownsProlongation\");\n");
    }
    writer.write("   // kernels::aderdg::optimised::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() ); //TODO JMG, uncomment in Toolkit when kernel implemented \n"); //TODO JMG, uncomment when kernel implemented
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"faceUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"faceUnknownsRestriction\");\n");
    }
    writer.write("  // kernels::aderdg::optimised::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() ); //TODO JMG, uncomment in Toolkit when kernel implemented \n"); //TODO JMG, uncomment when kernel implemented
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"faceUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"volumeUnknownsProlongation\");\n");
    }
    writer.write("  // kernels::aderdg::optimised::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() ); //TODO JMG, uncomment in Toolkit when kernel implemented \n"); //TODO JMG, uncomment when kernel implemented
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"volumeUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("   _profiler->start(\"volumeUnknownsRestriction\");\n");
    }
    writer.write("  // kernels::aderdg::optimised::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() ); //TODO JMG, uncomment in Toolkit when kernel implemented \n"); //TODO JMG, uncomment when kernel implemented
    if (_enableProfiler) {
      writer.write("   _profiler->stop(\"volumeUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS, double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {\n");
    if (_enableProfiler) {
        writer.write("  _profiler->start(\"boundaryConditions\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write(" // kernels::aderdg::generic::c" 
            + "::boundaryConditions"
            + "( *this, fluxOut, stateOut, fluxIn, stateIn, cellCentre, cellSize, t, dt, faceIndex, normalNonZero ); //TODO JMG, uncomment in Toolkit when kernel implemented \n"); //TODO JMG, uncomment when kernel implemented
    if (_enableProfiler) {
        writer.write("  _profiler->stop(\"boundaryConditions\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  /*
  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(
        solverName, writer, projectName, _numberOfUnknowns, _numberOfParameters, _order, _hasConstants);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    if (_dimensions == 2) {
      writer.write("void " + projectName + "::" + solverName
          + "::flux(const double* const Q, double* f, double* g) {\n");
    } else {
      writer.write("void " + projectName + "::" + solverName
          + "::flux(const double* const Q, double* f, double* g, double* h) {\n");
    }
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // f\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  f[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("  // g\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  g[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    if (_dimensions == 3) {
      writer.write("  // h\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  h[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a PDF.f90.\n");
  }
  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }
  */ //TODO JMG
  
    @Override
  public final void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(solverName, writer, projectName,
        _numberOfUnknowns, _numberOfParameters, _order, _hasConstants);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

      // flux
      writer.write("void " + projectName + "::" + solverName
            + "::flux(const double* const Q, double** F) {\n");
      writer.write("  // Dimensions             = " + _dimensions + "\n");
      writer.write(
          "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
              + " (#unknowns + #parameters)\n\n");
      writer.write("  double* f = F[0];\n");
      writer.write("  double* g = F[1];\n");
      if (_dimensions == 3) {
        writer.write("  double* h = F[2];\n");
      }
      writer.write("\n");
      writer.write("  // @todo Please implement\n");
      writer.write("  // f\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  f[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
      writer.write("  // g\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  g[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
      if (_dimensions == 3) {
        writer.write("  // h\n");
        writer.write("  // @todo Please implement\n");
        for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
          writer.write("  h[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
        }
      }
      writer.write("}\n");

    writer.write("\n\n\n");


      // source
      writer.write("void " + projectName + "::" + solverName + "::source(const double* const Q, double* S) {\n");
      writer.write("  // Number of variables = " + _numberOfUnknowns + " + " +  _numberOfParameters + "\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  S[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");
    writer.write("\n\n\n");
      
    // boundary conditions
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
            "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n\n");
    writer.write("\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  // fluxOut\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  fluxOut[" + String.format("%" + digits + "d", i) + "] = fluxIn[" + String.format("%" + digits + "d", i) + "];\n");
    }
    writer.write("  // stateOut\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  stateOut[" + String.format("%" + digits + "d", i) + "] = stateIn[" + String.format("%" + digits + "d", i) + "];\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    

      // eigenvalues
      writer.write("void " + projectName + "::" + solverName
          + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
      writer.write("  // Dimensions             = " + _dimensions + "\n");
      writer.write(
          "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
              + " (#unknowns + #parameters)\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
      writer.write("}\n");

    writer.write("\n\n\n");
    
    //initial conditions
    writer.write("bool " + projectName + "::" + solverName
        + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t, double dt) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return false;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
        "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    // refinement control
    writer.write("exahype::solvers::Solver::RefinementControl " + projectName + "::" + solverName
        + "::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return exahype::solvers::Solver::RefinementControl::Keep;\n");
    writer.write("}\n");
    writer.write("\n\n\n");

      // ncp

      writer.write("void " + projectName + "::" + solverName
          + "::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {\n");
      writer.write("  // Dimensions             = " + _dimensions + "\n");
      writer.write("  // Number of variables    = "
          + Integer.toString(_numberOfUnknowns + _numberOfParameters)
          + " (#unknowns + #parameters)\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _dimensions * (_numberOfUnknowns + _numberOfParameters); i++) {
          writer.write("  BgradQ[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");

    writer.write("\n\n\n");


      // matrixb
      writer.write("void " + projectName + "::" + solverName + "::matrixb(const double* const Q, const int normalNonZero, double* Bn) {\n");
      writer.write("  // Number of variables    = "
          + Integer.toString(_numberOfUnknowns + _numberOfParameters)
          + " (#unknowns + #parameters)\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < (_numberOfUnknowns + _numberOfParameters) * (_numberOfUnknowns + _numberOfParameters); i++) {
        writer.write("Bn[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");

    writer.write("\n\n\n");
  }
    
    @Override
    public boolean supportsVariables() {
      return false;
    }
}
