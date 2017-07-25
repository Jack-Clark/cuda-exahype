package eu.exahype.variables;

import java.io.BufferedWriter;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import eu.exahype.IOUtils;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.ANamingScheme;
import eu.exahype.node.ANamingSchemes;
import eu.exahype.node.AVariables;
import eu.exahype.node.AWithNameVariable;
import eu.exahype.node.AWithoutNameVariable;
import eu.exahype.node.PNamingScheme;
import eu.exahype.node.PSolver;
import eu.exahype.node.PVariable;

public class Variables {
  Map<String,Integer>             _variablesMap;
  int                             _numberOfVariables;
  Map<String,Integer>             _parametersMap;
  int                             _numberOfParameters;
  Map<String,Map<String,Integer>> _namingSchemesMap;
  
  public Map<String, Integer> getVariablesMap() {
    return _variablesMap;
  }

  public int getNumberOfVariables() {
    return _numberOfVariables;
  }

  public Map<String, Integer> getParametersMap() {
    return _parametersMap;
  }

  public int getNumberOfParameters() {
    return _numberOfParameters;
  }
  
  public Set<String> getNamingSchemeNames() {
    return _namingSchemesMap.keySet();
  }

  private static List<PVariable> getVariablesAsList(PSolver node) {
    List<PVariable> variablesAsList = null;
    if (node instanceof AAderdgSolver) {
      variablesAsList = ((AVariables)((AAderdgSolver) node).getVariables()).getVariable();
    } else if (node instanceof ALimitingAderdgSolver) {
      variablesAsList = ((AVariables)((ALimitingAderdgSolver) node).getVariables()).getVariable();
    } else if (node instanceof AFiniteVolumesSolver) {
      variablesAsList = ((AVariables)((AFiniteVolumesSolver) node).getVariables()).getVariable();
    } else {
      System.err.println("eu/exahype/variables/Variables.getVariablesAsList(PSolver)\t[ERROR]: Do not support solver type '"+node.getClass().toString()+"'.");
      System.exit(1);
    } 
    
    assert variablesAsList!=null;
    return variablesAsList;
  }
  
  /**
   * Returns an empty map if no parameters have been found.
   * Thus no check for null is necessary.
   */
  private static List<PVariable> getParametersAsList(PSolver node) {
    List<PVariable> parametersAsList = new LinkedList<PVariable>();
    
    if (node instanceof AAderdgSolver) {
      if ( ((AAderdgSolver) node).getParameters() != null) {
        parametersAsList.addAll(((AVariables)((AAderdgSolver) node).getParameters()).getVariable());
      }
    } else if (node instanceof ALimitingAderdgSolver) {
      if ( ((ALimitingAderdgSolver) node).getParameters() != null) {
        parametersAsList.addAll(((AVariables)((ALimitingAderdgSolver) node).getParameters()).getVariable());
      }
    } else if (node instanceof AFiniteVolumesSolver) {
      if ( ((AFiniteVolumesSolver) node).getParameters() != null) {
        parametersAsList.addAll(((AVariables)((AFiniteVolumesSolver) node).getParameters()).getVariable());
      }
    } else {
      System.out.println("eu/exahype/variables/Variables.getParametersAsList(PSolver)\t[WARNING]: Do not support solver type '"+node.getClass().toString()+"'.");
    }
    return parametersAsList;
  }
  
  private static Map<String,List<PVariable>> getNamingSchemesAsList(PSolver node) {
    List<PNamingScheme> namingSchemesAsList = null;
    
    if (node instanceof AAderdgSolver) {
      if ( ((AAderdgSolver) node).getNamingSchemes()!= null) {
        namingSchemesAsList = ((ANamingSchemes)(((AAderdgSolver) node).getNamingSchemes())).getNamingScheme();
      }
    } else if (node instanceof ALimitingAderdgSolver) {
      if ( ((ALimitingAderdgSolver) node).getNamingSchemes()!= null) {
        namingSchemesAsList = ((ANamingSchemes)(((ALimitingAderdgSolver) node).getNamingSchemes())).getNamingScheme();
      }
    } else if (node instanceof AFiniteVolumesSolver) {
      if ( ((AFiniteVolumesSolver) node).getNamingSchemes() != null) {
        namingSchemesAsList = ((ANamingSchemes)(((AFiniteVolumesSolver) node).getNamingSchemes())).getNamingScheme();
      }
    } else {
      System.out.println("eu/exahype/variables/Variables.getNamingSchemes(PSolver)\t[WARNING]: Do not support solver type '"+node.getClass().toString()+"'.");
    }
    
    Map<String,List<PVariable>> namingSchemes = new LinkedHashMap<String,List<PVariable>>(); // We want to preserve the ordering.
    if (namingSchemesAsList!=null) {
      for (PNamingScheme pNamingScheme : namingSchemesAsList) {
        ANamingScheme aNamingScheme = (ANamingScheme) pNamingScheme;
        String identifier = aNamingScheme.getName().getText();
        List<PVariable> variables = ((AVariables) aNamingScheme.getScheme()).getVariable();
        
        namingSchemes.put(identifier, variables);
      }
    }
    return namingSchemes;
  }
  
  private static Map<String,Integer> parseVariables(List<PVariable> variablesAsList,String defaultName) {
    Map<String, Integer> map = new LinkedHashMap<String, Integer>();
    
    for (PVariable pVariable :  variablesAsList) {
      if (pVariable instanceof AWithNameVariable) {
        AWithNameVariable variable = (AWithNameVariable) pVariable;
        String name       = variable.getName().getText();
        int multiplicity  = Integer.parseInt(variable.getMultiplicity().getText());
        
        if (multiplicity>0) {
          map.put(name, multiplicity);
        }
      } else if (pVariable instanceof AWithoutNameVariable) {
        AWithoutNameVariable variable = (AWithoutNameVariable) pVariable;
        String name       = defaultName;
        int multiplicity  = Integer.parseInt(variable.getMultiplicity().getText());
        
        if (multiplicity > 0) {
          map.put(name, multiplicity);
        }
      } else {
        System.out.println("ERROR: I do not know how to handle variable type "+pVariable.getClass().toString()+"!");
        System.exit(1);
        return null;
      }
    }
    return map;
  }
  
  /**
   * @note We rely on a linked hash map here that does not change the order of the variables.
   */
  public static Map<String, Integer> getVariables(PSolver node,String defaultName) {    
    List<PVariable> variablesAsList = getVariablesAsList(node);
    if (variablesAsList!=null) {
      return parseVariables(variablesAsList,defaultName); 
    } else {
      System.out.println("ERROR: I do not know how to handle variables of solver type "+node.getClass().toString()+"!");
      System.exit(1);
      return null;
    }
  }
  
  /**
   * @note We rely on a linked hash map here that does not change the order of the variables.
   */
  public static Map<String, Integer> getParameters(PSolver node,String defaultName) {    
    List<PVariable> variablesAsList = getParametersAsList(node);
    return parseVariables(variablesAsList,defaultName);
  }

  
  /**
   * @note We rely on a linked hash map here that does not change the order of the parameters.
   */
  public static Map<String, Map<String,Integer>> getNamingSchemes(PSolver node) {
    Map<String,Map<String,Integer>> namingSchemes = new LinkedHashMap<String,Map<String,Integer>>(); 
    
    Map<String,List<PVariable>> namingSchemesAsList = getNamingSchemesAsList(node);
    for (String name : namingSchemesAsList.keySet()) {
      namingSchemes.put(name, parseVariables(namingSchemesAsList.get(name), "values"));
    }
    return namingSchemes;
  }
  
  public static int sumMultiplicities(Map<String,Integer> variables) {
    int sumOfMultiplicities = 0;
    if (variables!=null) {
      for (String key : variables.keySet()) {
        sumOfMultiplicities += variables.get(key);
      }
    }
    
    return sumOfMultiplicities;
  }
  
  public Variables(PSolver node) {
    _variablesMap       = getVariables(node,"Q");
    _parametersMap      = getParameters(node,"params");
    _namingSchemesMap   = getNamingSchemes(node);
    _numberOfVariables  = sumMultiplicities(_variablesMap);
    _numberOfParameters = sumMultiplicities(_parametersMap);
  }
  
  public Variables(Map<String, Integer> variablesMap, Map<String, Integer> parametersMap, Map<String,Map<String, Integer>> namingSchemesMap, int dimensions) {
    _variablesMap       = variablesMap;
    _parametersMap      = parametersMap;
    _namingSchemesMap   = namingSchemesMap;
    _numberOfVariables  = sumMultiplicities(variablesMap);
    _numberOfParameters = sumMultiplicities(parametersMap);
  }
  
  private String appendVariableGetters(String getters, String identifier, int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double rho() const;
      getters += indent + "double "+identifier+"() const { return _Q["+offset+"]; }\n\n";
    } else {
      // ex: double v(int index) const;
      getters += indent + "double "+identifier+"(int index) const {\n"
              +  indent + "  assertion(index >= 0 && index<"+multiplicity+");\n"
              +  indent + "  return _Q["+offset+"+index];\n"
              +  indent + "}\n\n";
      // ex: tarch::la::Vector<3,double> v() const;
      getters += indent + "tarch::la::Vector<"+multiplicity+",double> "+identifier+"() const {\n"
              +  indent + "  tarch::la::Vector<"+multiplicity+",double> values;\n";
      getters += indent + "  values=";
      for (int i=0; i<multiplicity; i++) {
        getters += "_Q"+"["+(offset+i)+"]"+( i<multiplicity-1 ? "," : ";\n" );
      }
      getters += indent + "  return values;\n"
              +  indent + "}\n\n";
    }
    return getters;
  }
  
  /**
   * Generate the getters for the variables (and parameters).
   */
  private String createVariablesGetters() {
    String getters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      getters  = appendVariableGetters(getters,identifier,multiplicity,offset);
      offset  += multiplicity;
    }
    if (_parametersMap!=null) {
      for (String identifier : _parametersMap.keySet()) {
        int multiplicity = _parametersMap.get(identifier);
        getters  = appendVariableGetters(getters,identifier,multiplicity,offset);
        offset  += multiplicity;
      }
    }
    
    return getters;
  }
  
  private String appendVariableSetters(String setters, String identifier, int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      setters += indent + "double& "+identifier+"() { return _Q["+offset+"]; }\n\n";
    } else {
      // ex: double& v(int index);
      setters += indent +"double& "+identifier+"(int index) { return _Q["+offset+"+index]; }\n\n";
      
      // ex: void v(const tarch::la::Vector<3,double>& values);
      setters += indent +"void "+identifier+"(const tarch::la::Vector<"+multiplicity+",double>& values) {\n";
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"*(_Q+"+(offset+i)+")=values["+i+"];\n";
      }
      setters += indent +"}\n\n";
      
      // ex: void v(double v0, double v1, double v2);
      setters += indent +"void "+identifier+"(";
      for (int i=0; i<multiplicity; i++) {
        setters += "double "+identifier+i+( i<multiplicity-1 ? "," : ") {\n" );
      }
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"*(_Q+"+(offset+i)+")="+identifier+i+";\n";
      }
      setters += indent + "}\n\n";
    }
    return setters;
  }
  
  /**
   * Generate the setters for the variables (and parameters).
   */
  private String createVariablesSetters() {
    String setters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      setters = appendVariableSetters(setters, identifier, multiplicity, offset);
      offset += multiplicity;
    }
    if (_parametersMap!=null) {
      for (String identifier : _parametersMap.keySet()) {
        int multiplicity = _parametersMap.get(identifier);
        setters = appendVariableSetters(setters, identifier, multiplicity, offset);
        offset += multiplicity;
      }
    }
    
    return setters;
  }
  
  /**
   * Generate the getters for the primitives.
   */
  private String createNamingSchemeGetters(String name) {
    String getters = "";
    
    int offset = 0;
    if (_namingSchemesMap!=null) {
      for (String identifier : _namingSchemesMap.get(name).keySet()) {
        int multiplicity = _namingSchemesMap.get(name).get(identifier);
        getters  = appendVariableGetters(getters,identifier,multiplicity,offset);
        offset  += multiplicity;
      }
    }
    
    return getters;
  }
  
  /**
   * Generate the setters for the primitives.
   */
  private String createNamingSchemeSetters(String name) {
    String setters = "";
    
    int offset = 0;
    if (_namingSchemesMap!=null) {
      for (String identifier : _namingSchemesMap.get(name).keySet()) {
        int multiplicity = _namingSchemesMap.get(name).get(identifier);
        setters = appendVariableSetters(setters, identifier, multiplicity, offset);
        offset += multiplicity;
      }
    }
    
    return setters;
  }
  
  /**
   * @note The current implementation assumes a column major memory layout
   * of the fluxes.
   */
  private String appendFluxesGetters(String getters, String identifier,int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double rho(int column) const;  
      getters += indent + "double "+identifier+"(int column) const {\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"];\n"
              +  indent + "}\n\n";
      // ex: tarch::la::Vector<DIMENSIONS,double> rho() const;
      getters += indent + "tarch::la::Vector<DIMENSIONS,double> "+identifier+"() const {\n"
              +  indent + "  #if DIMENSIONS==2\n" 
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
              for (int i=0; i<2; i++) {
                getters += "_F["+i+"]["+offset+"]" + ( i<2-1 ? "," : ");\n" );
              }
     getters  += indent + "  #elif DIMENSIONS==3\n" 
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
              for (int i=0; i<3; i++) {
                getters += "_F["+i+"]["+offset+"]" + ( i<3-1 ? "," : ");\n" );
              }
      getters += indent + "  #endif\n"
              +  indent + "  return values;\n"
              +  indent + "}\n\n";
      
    } else {
      // ex: double v(int row, int column) const;  
      getters += indent + "double "+identifier+"(int row, int column) const {\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"+row];\n"
              +  indent + "}\n\n";
      
      // tarch::la::Vector<3,double> v(int row) const;
      getters += indent + "tarch::la::Vector<DIMENSIONS,double> "+identifier+"(int row) const {\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n"
              +  indent + "  #if DIMENSIONS==2\n"
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
              for (int i=0; i<2; i++) {
                getters += "_F["+i+"]["+offset+"+row]" + ( i<2-1 ? "," : ");\n" );
              }
      getters += indent + "  #elif DIMENSIONS==3\n" 
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
              for (int i=0; i<3; i++) {
                getters += "_F["+i+"]["+offset+"+row]" + ( i<3-1 ? "," : ");\n" );
              }
      getters += indent + "  #endif\n"
              +  indent + "  return values;\n"
              +  indent + "}\n\n";
      
      // ex: tarch::la::Matrix<3,3,double> v() const;
      getters += indent + "tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double> "+identifier+"() const {\n"
              +  indent + "  tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double> values;\n"
              +  indent + "  #if DIMENSIONS==2\n";
      getters += indent + "  values = ";      
              for (int i=0; i<multiplicity; i++) {
                for (int j=0; j<2; j++) {
                  getters += "_F["+j+"]["+(offset+i)+"]" + ( j<2-1 ? "," : "" );
                }
                getters += ( i<multiplicity-1 ? ",\n"+indent + "           " : ";\n" );
              }
      getters += indent + "  #elif DIMENSIONS==3\n";
      getters += indent + "  values = ";      
              for (int i=0; i<multiplicity; i++) {
                for (int j=0; j<2; j++) {
                  getters += "_F["+j+"]["+(offset+i)+"]" + ( j<2-1 ? "," : "" );
                }
                getters += ( i<multiplicity-1 ? ",\n"+indent + "           " : ";\n" );
              }
      getters += indent + "  #endif\n"
              +  indent + "  return values;\n"
              +  indent + "}\n\n";
    }
    return getters;
  }
  
  private String createFluxesGetters() {
    String getters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      getters  = appendFluxesGetters(getters,identifier,multiplicity,offset);
      offset  += multiplicity;
    }
    
    return getters;
  }

  private String appendFluxesSetters(String setters, String identifier,
      int multiplicity, int offset) {
    String indent  = "    ";

    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double& v(int column);
      setters += indent + "double& "+identifier+"(int column) {\n" 
          +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
          +  indent + "  return _F[column]["+offset+"];\n" 
          +  indent + "}\n\n";

      // ex: void rho(tarch::la::Vector<DIMENSIONS,double>& values); 3D and 2D
      setters += indent +"void "+identifier+"(const tarch::la::Vector<DIMENSIONS,double>& values) {\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"]=values["+j+"];\n";
      }
      setters += indent + "  #if DIMENSIONS==3\n";
      setters += indent + "  _F[2]["+offset+"]=values[2];\n";
      setters += indent + "  #endif\n";
      setters += indent + "}\n";
      // ex: void rho(tarch::la::Vector<DIMENSIONS,double>& values); 2.5D
      setters += indent + "#if DIMENSIONS==2\n";
      setters += indent +"/** Setter for 2.5D calculations. Third vector element is ignored.*/\n";
      setters += indent +"void "+identifier+"(const tarch::la::Vector<3,double>& values) {\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"]=values["+j+"];\n";
      }
      setters += indent + "}\n";
      setters += indent + "#endif\n\n";

      // ex: void v(double v0, double v1, double v2);  3D and 2.5D
      setters += indent +"/** Setter for 3D and 2.5D calculations. Third argument is ignored for the latter.*/\n";
      setters += indent +"void "+identifier+"(";
      for (int j=0; j<3; j++) {
        setters += "double v"+j+( j<2 ? "," : ") {\n" );
      }
      for (int j=0; j<2; j++) {
        setters += indent +"  _F["+j+"]["+offset+"]=v"+j+";\n";
      }
      setters += indent + "  #if DIMENSIONS==3\n";
      setters += indent + "  _F[2]["+offset+"]=v2;\n";
      setters += indent + "  #endif\n";
      setters += indent + "}\n";
      // ex: void v(double v0, double v1);  2D
      setters += indent + "#if DIMENSIONS==2\n";
      setters += indent + "void "+identifier+"(";
      for (int j=0; j<2; j++) {
        setters += "double v"+j+( j<2-1 ? "," : ") {\n" );
      }
      for (int j=0; j<2; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"]=v"+j+";\n";
      }
      setters += indent + "}\n";
      setters += indent + "#endif\n\n";
    } else {
      // ex: double& v(int row, int column);
      setters += indent + "double& "+identifier+"(int row, int column) {\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"+row];\n"
              +  indent + "}\n\n";

      // ex: void v(int row, const tarch::la::Vector<3,double>& values); 2D and 3D
      setters += indent +"void "+identifier+"(int row, const tarch::la::Vector<DIMENSIONS,double>& values) {\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  _F["+j+"]["+offset+"+row]=values["+j+"];\n";
      }
      setters += indent + "  #if DIMENSIONS==2\n";
      setters += indent +"  _F[2]["+offset+"+row]=values[2];\n";
      setters += indent + "  #endif\n";
      setters += indent +"}\n";
      setters += indent + "#if DIMENSIONS==2\n";
      setters += indent +"/** Setter for 2.5D calculations. Third vector element is ignored.*/\n";
      setters += indent +"void "+identifier+"(int row, const tarch::la::Vector<3,double>& values) {\n"
          +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  _F["+j+"]["+offset+"+row]=values["+j+"];\n";
      }
      setters += indent +"}\n";
      setters += indent +"#endif\n\n";
      
      // ex: void v(const tarch::la::Matrix<3,DIMENSIONS,double>& values);  2D and 3D
      setters += indent +"void "+identifier+"(const tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double>& values) {\n";
      for (int j=0; j<2; j++) {
        for (int i=0; i<multiplicity; i++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=values("+i+","+j+");\n";
        }
      }
      setters += indent + "  #if DIMENSIONS==3\n";
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"_F[2]["+(offset+i)+"]=values("+i+",2);\n";
      }
      setters += indent + "  #endif\n";
      setters += indent +"}\n";
      // ex: void v(const tarch::la::Matrix<3,DIMENSIONS,double>& values);  2.5
      setters += indent + "#if DIMENSIONS==2\n";
      setters += indent +"/** Setter for 2.5D calculations. Third matrix column is ignored.*/\n";
      setters += indent +"void "+identifier+"(const tarch::la::Matrix<"+multiplicity+",3,double>& values) {\n";
      for (int j=0; j<2; j++) {
        for (int i=0; i<multiplicity; i++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=values("+i+","+j+");\n";
        }
      }
      setters += indent +"}\n";
      setters += indent +"#endif\n\n";
      
      // ex: void v(int row, double v0, double v1, double v2);
      setters += indent +"/** Setter for 3D and 2.5D calculations. Third argument is ignored for the latter.*/\n";
      setters += indent +"void "+identifier+"(int row, ";
      for (int j=0; j<3; j++) {
        setters += "double v"+j+( j<2 ? "," : ") {\n" );
      }
      setters += indent + "  assertion(row >= 0 && row<"+multiplicity+");\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  _F["+j+"]["+offset+"+row]=v"+j+";\n";
      }
      setters += indent +"  #if DIMENSIONS==3\n";
      setters += indent +"  _F[2]["+offset+"+row]=v2;\n";
      setters += indent +"  #endif\n";
      setters += indent +"}\n";
      setters += indent +"#if DIMENSIONS==2\n";
      setters += indent +"/** Setter for 2D calculations.*/\n";
      setters += indent +"void "+identifier+"(int row, ";
      for (int j=0; j<2; j++) {
        setters += "double v"+j+( j<1 ? "," : ") {\n" );
      }
      setters += indent + "  assertion(row >= 0 && row<"+multiplicity+");\n";
      for (int j=0; j<2; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"+row]=v"+j+";\n";
      }
      setters += indent +"}\n";
      setters += indent +"#endif\n\n";

      // ex: void v(double v00, double v01, double v02, double v10, double v11, double v12, double v20, double v21, double v22);  2.5D and 3D
      int argIndentLength = ( indent +"void "+identifier+"(" ).length();
      String argIndent = new String(new char[argIndentLength]).replace("\0", " ");
      setters += indent +"/** Setter for 3D and 2.5D calculations. Third column values are ignored for the latter.*/\n";
      setters += indent +"void "+identifier+"(";
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<3; j++) {
          setters += "double v"+i+j+( j<2 ? ", " : "" );
        }
        setters += ( i<multiplicity-1 ? ",\n"+argIndent : ") {\n" );
      }
      for (int j=0; j<2; j++) {
        for (int i=0; i<multiplicity; i++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=v"+i+j+";\n";
        }
      }
      setters += indent +"  #if DIMENSIONS==3\n";
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"_F[2]["+(offset+i)+"]=v"+i+"2;\n";
      }
      setters += indent +"  #endif\n";
      setters += indent +"}\n";
      // ex: void v(double v00, double v01, double v02, double v10, double v11, double v12, double v20, double v21, double v22);  2D
      setters += indent +"#if DIMENSIONS==2\n";
      setters += indent +"void "+identifier+"(";
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<2; j++) {
          setters += "double v"+i+j+( j<1 ? ", " : "" );
        }
        setters += ( i<multiplicity-1 ? ",\n"+argIndent : ") {\n" );
      }
      for (int j=0; j<2; j++) {
        for (int i=0; i<multiplicity; i++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=v"+i+j+";\n";
        }
      }
      setters += indent +"}\n";
      setters += indent +"#endif\n\n";
    }
    return setters;
  }
  
  private String createFluxesSetters() {
    String setters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      setters = appendFluxesSetters(setters, identifier, multiplicity, offset);
      offset += multiplicity;
    }
    
    return setters;
  }

  public void writeHeader(BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/variables/templates/VariablesHeader.template");
    
    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}",  solverName);
    
    content = content.replaceAll("\\{\\{NumberOfVariables\\}\\}",  String.valueOf(_numberOfVariables));
    content = content.replaceAll("\\{\\{NumberOfParameters\\}\\}", String.valueOf(_numberOfParameters));
    
    String variablesGetters = createVariablesGetters();
    String variablesSetters = createVariablesSetters();
    content = content.replaceAll("\\{\\{VariablesGetters\\}\\}", variablesGetters);
    content = content.replaceAll("\\{\\{VariablesSetters\\}\\}", variablesSetters);
    
    String fluxesGetters = createFluxesGetters();
    String fluxesSetters = createFluxesSetters();
    content = content.replaceAll("\\{\\{FluxesGetters\\}\\}", fluxesGetters);
    content = content.replaceAll("\\{\\{FluxesSetters\\}\\}", fluxesSetters);
    
    // now generate the naming schemes
    String namingSchemesTemplate = IOUtils.convertRessourceContentToString(
        "eu/exahype/variables/templates/NamingScheme.template");
    
    String namingSchemes = "";
    for (String name : _namingSchemesMap.keySet()) {
      int multiplicities = sumMultiplicities(_namingSchemesMap.get(name));
      
      namingSchemes += namingSchemesTemplate+"\n\n";
      namingSchemes = namingSchemes.replaceAll("\\{\\{Project\\}\\}", projectName);
      namingSchemes = namingSchemes.replaceAll("\\{\\{Solver\\}\\}",  solverName);
      namingSchemes = namingSchemes.replaceAll("\\{\\{Name\\}\\}",  name.substring(0,1).toUpperCase()+name.substring(1));
      namingSchemes = namingSchemes.replaceAll("\\{\\{Size\\}\\}",  String.valueOf(multiplicities));
      
      String getters = createNamingSchemeGetters(name).replaceAll("_Q","_data");
      String setters = createNamingSchemeSetters(name).replaceAll("_Q","_data");
      namingSchemes = namingSchemes.replaceAll("\\{\\{Getters\\}\\}", getters);
      namingSchemes = namingSchemes.replaceAll("\\{\\{Setters\\}\\}", setters);
    }
    
    content = content.replaceAll("\\{\\{NamingSchemes\\}\\}", namingSchemes);
    writer.write(content);
  }
}
