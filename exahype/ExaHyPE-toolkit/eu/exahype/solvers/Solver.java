package eu.exahype.solvers;

public interface Solver {
 
  /**
   * @return true if the solver supports generation of Variables classes.
   */
  public boolean supportsVariables();
  
  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException;
  
  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException;
  
  /**
   * @deprecated This method will be removed as soon as an abstract base class
   * is generated for the optimised solvers.
   */
  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException;
  
  public void writeAbstractHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException;
  
  public void writeAbstractImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException;
}
