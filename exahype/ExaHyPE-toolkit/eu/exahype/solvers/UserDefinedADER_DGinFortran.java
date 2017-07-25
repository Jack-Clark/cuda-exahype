package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class UserDefinedADER_DGinFortran implements Solver {
  public static final String Identifier = UserDefinedADER_DGinC.Identifier;

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
    // @todo
    System.err.println("not implemented yet\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    // @todo Implement
    System.err.println("not implemented yet\n");
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
