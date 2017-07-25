package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class OptimisedFluxesLinearADER_DGinC implements Solver {
  public static final String Identifier = "optimised::fluxes::linear";

  private int _dimensions;
  private int _numberOfUnknowns;
  private int _numberOfParameters;
  private int _order;
  private String _microarchitecture;
  private String _pathToLibxsmm;
  private boolean _enableProfiler;
  private boolean _hasConstants;

  public OptimisedFluxesLinearADER_DGinC(int dimensions, int numberOfUnknowns, int numberOfParameters, int order,
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
    Helpers.invokeCodeGenerator(solverName, _numberOfUnknowns, _numberOfParameters, _order, true, _dimensions,
        _microarchitecture, _pathToLibxsmm);

    writer.write("  private:\n");
    if (_dimensions == 2) {
      writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
    } else {
      writer.write(
          "    static void flux(const double* const Q, double* f, double* g, double* h);\n");
    }
    writer.write(
        "    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write(
        "    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");

    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.invokeCodeGenerator(solverName, _numberOfUnknowns, _numberOfParameters, _order, true, _dimensions,
        _microarchitecture, _pathToLibxsmm);

    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/optimised/Kernels.h\"\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {\n");
    // Cauchy-Kowalewski
    if (_enableProfiler) {
      writer.write("    _profiler->start(\"spaceTimePredictor\");\n");
    }
    writer.write("   kernels::aderdg::optimised::predictor( lQhi, lFhi, lQi, lFi );\n");
    writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, lQhi, lFhi );\n");
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
    // @todo boundaryConditions are missing
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
  }

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
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters) + " (#unknowns + #parameters)\n");
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
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters) + " (#unknowns + #parameters)\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters) + " (#unknowns + #parameters)\n");
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
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
