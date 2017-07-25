package eu.exahype.solvers;

import java.util.Set;

public class SolverFactory {
  private int _dimensions;
  private boolean _enableProfiler;
  private String _microarchitecture;
  private String _pathToLibxsmm;

  public SolverFactory(
      int dimensions,
      boolean enableProfiler,
      String microarchitecture,
      String pathToLibxsmm) {
    _dimensions = dimensions;
    _enableProfiler = enableProfiler;
    _microarchitecture = microarchitecture;
    _pathToLibxsmm = pathToLibxsmm;
  }
  
  public Solver createADERDGSolver(String kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int order,boolean hasConstants) {
    String generalKernel = kernel.substring(0, kernel.lastIndexOf("::"));
    boolean isLinear     = kernel.substring(kernel.lastIndexOf("::")).equalsIgnoreCase("::linear");
    
    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinFortran.Identifier )) {
      return new eu.exahype.solvers.UserDefinedADER_DGinFortran();
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinC.Identifier )) {
      return new eu.exahype.solvers.UserDefinedADER_DGinC(numberOfVariables,
          numberOfParameters, order, hasConstants, _enableProfiler);
    }
    else if (generalKernel.equals( eu.exahype.solvers.GenericADERDG.Identifier )) {
      return new eu.exahype.solvers.GenericADERDG(_dimensions,
          numberOfVariables, numberOfParameters, namingSchemeNames, order, _enableProfiler, hasConstants, isLinear, isFortran );
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )) {
      return new eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC(_dimensions,
          numberOfVariables, numberOfParameters, order, _microarchitecture, _pathToLibxsmm,
          _enableProfiler, hasConstants);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedADERDG.Identifier )) {
      return new eu.exahype.solvers.OptimisedADERDG(_dimensions,
          numberOfVariables, numberOfParameters, order, _microarchitecture, _pathToLibxsmm,
          _enableProfiler, hasConstants, false);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.KernelEuler2d.Identifier )) {
      return new eu.exahype.solvers.KernelEuler2d();
    }

    return null;
  }
  
  public Solver createFiniteVolumesSolver(String kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int patchSize,boolean hasConstants) {
    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinFortran.Identifier )) {
      return new eu.exahype.solvers.UserDefinedFiniteVolumesinFortran(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinC.Identifier )) {
      return new eu.exahype.solvers.UserDefinedFiniteVolumesinC(_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesGodunovInC.Identifier )) {
      return new eu.exahype.solvers.GenericFiniteVolumesGodunovInC(_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
    }

    return null;
  }
}
