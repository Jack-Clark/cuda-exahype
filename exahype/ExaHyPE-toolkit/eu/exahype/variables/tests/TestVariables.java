package eu.exahype.variables.tests;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.LinkedHashMap;
import java.util.Map;

import eu.exahype.variables.Variables;

public class TestVariables {
  
  public static void testVariables() {
    int dimensions=2;
    
    Map<String,Integer> variablesMap = new LinkedHashMap<String,Integer>();
    variablesMap.put("rho", 1);
    variablesMap.put("j", 3);
    variablesMap.put("E", 1);
    
    Map<String,Map<String,Integer>> namingSchemesMap = new LinkedHashMap<String,Map<String,Integer>>();
    Map<String,Integer> primitivesMap = new LinkedHashMap<String,Integer>();
    primitivesMap.put("rho", 1);
    primitivesMap.put("u", 3);
    primitivesMap.put("E", 1);
    namingSchemesMap.put("primitives", primitivesMap);
    
    Map<String,Integer> tools = new LinkedHashMap<String,Integer>();
    tools.put("hammer", 1);
    tools.put("nails", 3);
    tools.put("screwdriver", 10);
    namingSchemesMap.put("tools", tools);
    
    Map<String,Integer> parametersMap = new LinkedHashMap<String,Integer>();
    parametersMap.put("matScalar", 1);
    parametersMap.put("matVector", 3);
    
    Variables variables = new Variables(variablesMap,parametersMap,namingSchemesMap,dimensions);
    
    BufferedWriter bufferedWriter = new BufferedWriter(new OutputStreamWriter(System.out));
    
    try {
      System.out.println("Content of header file: ");
      variables.writeHeader(bufferedWriter, "MySolver", "MyProject");
      bufferedWriter.flush();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  public static void main(String[] args) {
    testVariables();
  }
}
