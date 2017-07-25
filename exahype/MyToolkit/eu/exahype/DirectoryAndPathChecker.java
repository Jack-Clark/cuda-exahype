package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;
import eu.exahype.node.APaths;

public class DirectoryAndPathChecker extends DepthFirstAdapter {
  public Boolean valid = true;

  protected java.io.File peanoKernelPath;
  protected java.io.File peanoToolboxPath;
  protected java.io.File exahypePath;
  protected java.io.File outputDirectory;
  protected java.io.File sharedMemoryOraclesPath;

  @Override
  public void inAPaths(APaths node) {
    peanoKernelPath  = new java.io.File(node.getPeanoKernelPath().getText());
    peanoToolboxPath = node.getPeanoKernelPath()==null ? peanoKernelPath : new java.io.File(node.getPeanoKernelPath().getText());
    exahypePath = new java.io.File(node.getExahypePath().getText());
    outputDirectory = new java.io.File(node.getOutputDirectory().getText());

    System.out.print("Peano kernel path: " + peanoKernelPath.getAbsolutePath());
    if (peanoKernelPath.isDirectory()) {
      System.out.println(" ... ok");
    } else {
      System.out.println(" ... not found");
      valid = false;
    }

    System.out.print("Peano kernel path " + peanoKernelPath.getAbsolutePath() + " holds peano sources");
    if ((new java.io.File(peanoKernelPath.getAbsolutePath() + "/peano")).isDirectory()) {
        System.out.println(" ... ok");
      } else {
        System.out.println(" ... not found");
        valid = false;
      }

    System.out.print("Peano kernel path " + peanoKernelPath.getAbsolutePath() + " holds tarch sources");
    if ((new java.io.File(peanoKernelPath.getAbsolutePath() + "/tarch")).isDirectory()) {
        System.out.println(" ... ok");
      } else {
        System.out.println(" ... not found");
        valid = false;
      }

    System.out.print("Peano toolboxes path: " + peanoToolboxPath.getAbsolutePath());
    if (peanoToolboxPath.isDirectory()) {
      System.out.println(" ... ok");
    } else {
      System.out.println(" ... not found");
      valid = false;
    }

    System.out.print("Peano toolboxes path " + peanoToolboxPath.getAbsolutePath() + " holds multiscalelinkedcell sources");
    if ((new java.io.File(peanoToolboxPath.getAbsolutePath() + "/multiscalelinkedcell")).isDirectory()) {
        System.out.println(" ... ok");
      } else {
        System.out.println(" ... not found");
        valid = false;
      }

    System.out.print("Peano toolboxes path " + peanoToolboxPath.getAbsolutePath() + " holds sharedmemoryoracles sources");
    if ((new java.io.File(peanoToolboxPath.getAbsolutePath() + "/sharedmemoryoracles")).isDirectory()) {
        System.out.println(" ... ok");
      } else {
        System.out.println(" ... not found");
        valid = false;
      }

    System.out.print("Peano toolboxes path " + peanoToolboxPath.getAbsolutePath() + " holds mpibalancing sources");
    if ((new java.io.File(peanoToolboxPath.getAbsolutePath() + "/mpibalancing")).isDirectory()) {
        System.out.println(" ... ok");
      } else {
        System.out.println(" ... not found");
        valid = false;
      }

  };

  @Override
  public void outAProject(AProject node) {
    System.out.print("ExaHyPE path: " + exahypePath.getAbsolutePath());
    if (exahypePath.isDirectory()) {
      System.out.println(" ... ok");
    } else {
      System.out.println(" ... not found");
      valid = false;
    }

    System.out.print("output directory: " + outputDirectory.getAbsolutePath());
    if (outputDirectory.isDirectory()) {
      System.out.println(" ... does exist (will not be overwritten)");
    } else {
      boolean createdDirectories = outputDirectory.mkdirs();

      if (createdDirectories) {
        System.out.println(" ... created");
      } else {
        System.out.println(" ... not found and could not be created");
        valid = false;
      }
    }
  }
}
