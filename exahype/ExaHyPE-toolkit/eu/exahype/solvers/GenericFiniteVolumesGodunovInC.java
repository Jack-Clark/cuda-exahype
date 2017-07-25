package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Set;

import eu.exahype.IOUtils;

public class GenericFiniteVolumesGodunovInC implements Solver {
  public static final String Identifier = "generic::Godunov";

  private int _dimensions;
  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _patchSize;
  private Set<String> _namingSchemeNames;
//  private int _patchSize;
  private boolean _enableProfiler;
  private boolean _hasConstants;

  public GenericFiniteVolumesGodunovInC(int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames, int patchSize,
      boolean enableProfiler, boolean hasConstants) {
    _dimensions = dimensions;
    _numberOfVariables = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _namingSchemeNames = namingSchemeNames;
    _patchSize = patchSize;
    _enableProfiler = enableProfiler;
    _hasConstants = hasConstants;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericFiniteVolumesSolverHeader.template");

    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);

    String profilerInclude                     = "";
    String solverConstructorSignatureExtension = "";
    if (_enableProfiler) {
      profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
      solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler"; 
    }
    if (_hasConstants) {
      solverConstructorSignatureExtension += ", exahype::Parser::ParserView constants"; // TODO(Dominic): Why pass by value? 
    }


    content = content.replaceAll("\\{\\{NumberOfVariables\\}\\}", String.valueOf(_numberOfVariables));
    content = content.replaceAll("\\{\\{NumberOfParameters\\}\\}",String.valueOf( _numberOfParameters));
    content = content.replaceAll("\\{\\{Dimensions\\}\\}",String.valueOf( _dimensions));
    //content = content.replaceAll("\\{\\{Order\\}\\}", String.valueOf(_order)); // Goudonov is 2nd order or so. Should probably tell here.

    content = content.replaceAll("\\{\\{ProfilerInclude\\}\\}",profilerInclude);
    content = content.replaceAll("\\{\\{SolverConstructorSignatureExtension\\}\\}", solverConstructorSignatureExtension);

    writer.write(content);
  }
  
  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericFiniteVolumesSolverGodunovInCUserCode.template");
    
    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);
    
    content = content.replaceAll("\\{\\{Elements\\}\\}",  String.valueOf( _numberOfParameters+_numberOfVariables));
    content = content.replaceAll("\\{\\{Dimensions\\}\\}",String.valueOf(_dimensions));
    
    //    String SolverInitSignatureExtension = "";
    String SolverInitSignatureExtension = "";
    if (_hasConstants) {
      SolverInitSignatureExtension = ", exahype::Parser::ParserView& constants";
    }
    content = content.replaceAll("\\{\\{SolverInitSignatureExtension\\}\\}", SolverInitSignatureExtension);
    
    // 
    int digits = String.valueOf(_numberOfVariables + _numberOfParameters).length();

    String adjustedSolutionValues = "  // State variables:\n";
    for (int i = 0; i < _numberOfVariables; i++) {
      adjustedSolutionValues += "  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) adjustedSolutionValues += "\n";
    }
    if (_numberOfParameters>0) {
      adjustedSolutionValues += "  // Material parameters:\n";
      for (int i = 0; i < _numberOfParameters; i++) {
        adjustedSolutionValues += "  Q[" + String.format("%" + digits + "d", _numberOfVariables+i) + "] = 0.0;";
        if (i<_numberOfParameters-1) adjustedSolutionValues += "\n";
      }
    }
    String SolverInitCallExtension             = "";
    if (_hasConstants) {
       SolverInitCallExtension = ", constants";
    }
    
    content = content.replaceAll("\\{\\{SolverInitCallExtension\\}\\}",SolverInitCallExtension);

    String eigenvalues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      eigenvalues += "  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) eigenvalues += "\n";
    }

    String flux = "";
    for (int d=0; d<_dimensions; ++d) {
      for (int i = 0; i < _numberOfVariables; i++) {
        flux += "  F["+d+"][" + String.format("%" + digits + "d", i) + "] = 0.0;";
        if (i<_numberOfVariables-1) flux += "\n";
      }
      if (d<_dimensions-1) {
        flux += "\n\n";    
      }
    }

    String source = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      source += "  S[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) source += "\n";
    }
    
    String boundaryValues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      boundaryValues += "  stateOutside[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) boundaryValues += "\n";
    }
    
    content = content.replaceAll("\\{\\{AdjustedSolutionValues\\}\\}",adjustedSolutionValues);
    content = content.replaceAll("\\{\\{Eigenvalues\\}\\}",eigenvalues);
    content = content.replaceAll("\\{\\{Flux\\}\\}",flux);
    content = content.replaceAll("\\{\\{Source\\}\\}",source);
    content = content.replaceAll("\\{\\{BoundaryValues\\}\\}",boundaryValues);
    
    writer.write(content);
  }

  /**
   * @deprecated This will be removed after the optimised kernel code generation
   * supports the abstract solvers.
   */
  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/GenericFiniteVolumesSolverGodunovInCGeneratedCode.template");
    
      content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
      content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);
      //
      String profilerInclude                     = "";
      String solverConstructorSignatureExtension = "";
      String solverConstructorArgumentExtension  = "";
      if (_enableProfiler) {
          profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
          solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler";
          solverConstructorArgumentExtension  += ", std::move(profiler)";
      }
      if (_hasConstants) {
          solverConstructorSignatureExtension += ", exahype::Parser::ParserView constants"; // TODO(Dominic): Why pass by value? 
      }
      content = content.replaceAll("\\{\\{ProfilerInclude\\}\\}",profilerInclude);
      content = content.replaceAll("\\{\\{SolverConstructorSignatureExtension\\}\\}", solverConstructorSignatureExtension);
      content = content.replaceAll("\\{\\{SolverConstructorArgumentExtension\\}\\}", solverConstructorArgumentExtension);
      //
      content = content.replaceAll("\\{\\{NumberOfVariables\\}\\}", String.valueOf(_numberOfVariables));
      content = content.replaceAll("\\{\\{NumberOfParameters\\}\\}",String.valueOf( _numberOfParameters));
      
      writer.write(content);
  }
  
  @Override
  public void writeAbstractHeader(BufferedWriter writer, String solverName,
      String projectName) throws IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/AbstractGenericFiniteVolumesSolverHeader.template");

    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);

    String profilerInclude                     = "";
    String solverConstructorSignatureExtension = "";
    if (_enableProfiler) {
      profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
      solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler"; 
    }

    content = content.replaceAll("\\{\\{ProfilerInclude\\}\\}",profilerInclude);
    content = content.replaceAll("\\{\\{SolverConstructorSignatureExtension\\}\\}", solverConstructorSignatureExtension);
    
    content = content.replaceAll("\\{\\{NumberOfVariables\\}\\}", String.valueOf(_numberOfVariables));
    content = content.replaceAll("\\{\\{NumberOfParameters\\}\\}",String.valueOf( _numberOfParameters));
    content = content.replaceAll("\\{\\{Dimensions\\}\\}",String.valueOf( _dimensions));
    content = content.replaceAll("\\{\\{PatchSize\\}\\}", String.valueOf(_patchSize));
    
    String namingSchemes = "";
    for (String name : _namingSchemeNames) {
      namingSchemes += "    " + "class "+name.substring(0, 1).toUpperCase() + name.substring(1) + ";\n";
    }
    content = content.replaceAll("\\{\\{NamingSchemes\\}\\}", namingSchemes);

    writer.write(content);
  }
  
  @Override
  public void writeAbstractImplementation(BufferedWriter writer,
      String solverName, String projectName) throws IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/solvers/templates/AbstractGenericFiniteVolumesSolverGodunovInCImplementation.template");

    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);
    //
    String profilerInclude                     = "";
    String solverConstructorSignatureExtension = "";
    String solverConstructorArgumentExtension  = "";
    if (_enableProfiler) {
      profilerInclude                        = "#include \"exahype/profilers/Profiler.h\"";
      solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler";
      solverConstructorArgumentExtension  += ", std::move(profiler)";
    }
    if (_hasConstants) {
      solverConstructorSignatureExtension += ", exahype::Parser::ParserView constants"; // TODO(Dominic): Why pass by value? 
    }
    content = content.replaceAll("\\{\\{ProfilerInclude\\}\\}",profilerInclude);
    content = content.replaceAll("\\{\\{SolverConstructorSignatureExtension\\}\\}", solverConstructorSignatureExtension);
    content = content.replaceAll("\\{\\{SolverConstructorArgumentExtension\\}\\}", solverConstructorArgumentExtension);
    //
    content = content.replaceAll("\\{\\{NumberOfVariables\\}\\}", String.valueOf(_numberOfVariables));
    content = content.replaceAll("\\{\\{NumberOfParameters\\}\\}",String.valueOf( _numberOfParameters));

    // TODO(Dominic): Add profilers
    
    writer.write(content); 
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
    return true;
  }
}
