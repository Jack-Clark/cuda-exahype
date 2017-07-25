package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class KernelEuler2d implements Solver {
  public static final String Identifier = "kernel::euler2d";

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo
    System.err.println("not implemented yet\n");
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
  
  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/optimised/Kernels.h\"\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor(double* lQhbnd, double* lFhbnd, double** tempSpaceTimeUnknowns, double** tempSpaceTimeFluxUnknowns, double*  tempUnknowns,double*  tempFluxUnknowns, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt) {\n");
    writer.write("   kernels::aderdg::optimised::picardLoop<flux>( tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0], luh, dx, dt );\n");
    writer.write("   kernels::aderdg::optimised::predictor( tempUnknowns, tempFluxUnknowns, tempSpaceTimeUnknowns[0], tempSpaceTimeFluxUnknowns[0] );\n");
    writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, tempUnknowns, tempFluxUnknowns );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate( double* luh, const double* const lduh, const double dt ) {\n");
    writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral( double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
    writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral( double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
    writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver( double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex ) {\n");
    writer.write(
        "   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    // @todo boundary conditions are missing
    //
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize( const double* const luh, double* tempEigenvalues, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
    writer.write(
        "   return kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n"); // todo do not pass the temp array yet.
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment( double*  luh, const tarch::la::Vector<DIMENSIONS,double>&   center, const tarch::la::Vector<DIMENSIONS,double>&   dx, double  t, double  dt ) {\n");
    writer.write(
        "   kernels::aderdg::generic::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
  }
  
  @Override
  public boolean supportsVariables() {
    return false;
  }
}
