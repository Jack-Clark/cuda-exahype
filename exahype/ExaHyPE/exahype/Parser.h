/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifndef EXAHYPE_PARSER
#define EXAHYPE_PARSER

namespace exahype {
class Parser;
}

#include <iostream>

#include <map>
#include <vector>

#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "exahype/solvers/Solver.h"

/**
 * ExaHyPE command line parser
 *
 * A simple parser that creates a linear token stream
 * from a given ExaHyPE specification file.
 *
 * The parser can deal with C and doxygen style comments as
 * long as not two multi-line comment blocks are opened/closed
 * in the same line.
 *
 * @author Tobias Weinzierl, Dominic Etienne Charrier
 */
class exahype::Parser {
 public:
  /**
   * View on the parser
   *
   * An instance of this class is a parser view. While the parser sees the
   * whole specification file, a view only 'sees' the fragment that is
   * specific to one solver. As such, we do pass it to solvers (that hold
   * constants) and then allow these solvers to read their data from the
   * config file.
   *
   * From the user's point of view, this class provides access to key-value
   * pairs. If you have an instruction alike
   * <pre>
  constants         = {rho:0.4567,gamma:-4,alpha:4.04e-5,file:output}
     </pre>
   * in your
   *
   * @author Tobias Weinzierl
   */
  class ParserView {
   private:
    Parser& _parser;
    int _solverNumberInSpecificationFile;

    /**
     * @return Value for given key. Returns the empty string if there is no
     *         value. Returns std::npos if the key does not exist.
     */
    std::string getValue(const std::string selector,
                         const std::string& key) const;

   public:
    ParserView(Parser& parser, int solverNumberInSpecificationFile);

    /**
     * You may use keys without a value. This operation allows you to check
     * whether there are such keys. Furthermore, you might use this guy as
     * a preamble to the other getters.
     */
    bool hasKey(const std::string& key) const;

    /**
     * Please ensure that isValueValidXXX holds before you invoke this
     * operation.
     */
    bool getValueAsBool(const std::string& key) const;

    /**
     * Please ensure that isValueValidXXX holds before you invoke this
     * operation.
     */
    int getValueAsInt(const std::string& key) const;

    /**
     * Please ensure that isValueValidXXX holds before you invoke this
     * operation.
     */
    double getValueAsDouble(const std::string& key) const;

    /**
     * Please ensure that isValueValidXXX holds before you invoke this
     * operation.
     */
    std::string getValueAsString(const std::string& key) const;

    bool isValueValidBool(const std::string& key) const;
    bool isValueValidInt(const std::string& key) const;
    bool isValueValidDouble(const std::string& key) const;
    bool isValueValidString(const std::string& key) const;
  };

 private:
  static tarch::logging::Log _log;

  static const std::string   _noTokenFound;

  std::vector<std::string> _tokenStream;

  /*
   * Helper map for converting strings to types.
   */
  std::map<std::string, exahype::solvers::Solver::Type> _identifier2Type;

  /*
   * Helper map for converting strings to types.
   */
  std::map<std::string, exahype::solvers::Solver::TimeStepping>
      _identifier2TimeStepping;

  /**
   * Has to be static. If it is not static, then we can't modify it inside
   * const functions, i.e. all getters have to become non-const. This would
   * be reasonable but then in turn enforce all operations accepting parsers
   * to accept them as non-const.
   */
  static bool _interpretationErrorOccured;

  /**
   * \return "notoken" if not found.
   */
  std::string getTokenAfter(std::string token,
                            int additionalTokensToSkip = 0) const;
  std::string getTokenAfter(std::string token0, std::string token1,
                            int additionalTokensToSkip = 0) const;
  std::string getTokenAfter(std::string token0, int occurance0,
                            int additionalTokensToSkip) const;
  std::string getTokenAfter(std::string token0, int occurance0,
                            std::string token1, int occurance1,
                            int additionalTokensToSkip = 0) const;

 public:
  /**
   * Property strings in ExaHyPE are string alike "{all,left=0.5,Q4}". This
   * operation returns the value of a property, i.e. if you invoke
   * getvalueFromProperyString( "left" ), you obtain 0.5 in the example
   * above. The routine returns nan if now entry is found or the entry's
   * value  is not a valid floating point number.
   */
  static double getValueFromPropertyString(const std::string& parameterString,
                                           const std::string& key);

  Parser();
  virtual ~Parser() {}

  // Disallow copy and assignment
  Parser(const Parser& other) = delete;
  Parser& operator=(const Parser& other) = delete;

  enum class MulticoreOracleType {
    Dummy,
    Autotuning,
    GrainSizeSampling
    // evtl. spaeter mal InvadeSHM
  };

  enum class MPILoadBalancingType { Static };

  /**
   * <h2>Limitations</h2>
   * Requires a proper implementation of the
   * C++11 regex routines (GCC>=4.9.0).
   *
   * Cannot deal with multiple openings/closings of comments
   * in the same line.
   */
  void readFile(const std::string& filename);

  bool isValid() const;

  /**
   * \return How many threads is the code supposed to use?
   */
  int getNumberOfThreads() const;

  tarch::la::Vector<DIMENSIONS, double> getDomainSize() const;

  /**
   * @return Bounding box size. If the user has specified a non-cubical domain,
   *         then the bounding box still is cubical and all of its entries are
   *         the biggest dimension along one coordinate axis.
   */
  tarch::la::Vector<DIMENSIONS, double> getBoundingBoxSize() const;

  tarch::la::Vector<DIMENSIONS, double> getOffset() const;

  std::string getMulticorePropertiesFile() const;

  MulticoreOracleType getMulticoreOracleType() const;

  MPILoadBalancingType getMPILoadBalancingType() const;
  std::string getMPIConfiguration() const;
  int getMPIBufferSize() const;
  int getMPITimeOut() const;

  double getSimulationEndTime() const;

  /**
   * \return Indicates if the user has chosen the fused ADER-DG time stepping
   * variant.
   *
   * If the parser returns _noTokenFound, we may not issue an error as this is
   * an optional entry in the spec file.
   */
  bool getFuseAlgorithmicSteps() const;

  /**
   * \return Time step size underestimation factor for the fused ADER-DG time
   * stepping variant.
   */
  double getFuseAlgorithmicStepsFactor() const;

  double getTimestepBatchFactor() const;
  bool getSkipReductionInBatchedTimeSteps() const;

  double getDoubleCompressionFactor() const;
  bool   getSpawnDoubleCompressionAsBackgroundTask() const;

  /**
   * If we batch time steps, we can in principle switch off the boundary data
   * exchange, as ExaHyPE's data flow is realised through heaps. However, if we
   * turn off the boundary exchange, we enforce that no AMR is done in-between
   * domain boundaries.
   */
  bool getExchangeBoundaryDataInBatchedTimeSteps() const;

  /**
   * \return The type of a solver.
   */
  exahype::solvers::Solver::Type getType(int solverNumber) const;

  /**
   * \return The identifier of a solver.
   */
  std::string getIdentifier(int solverNumber) const;

  /**
   * \return The number of state vaParserriables of a solver.
   */
  int getVariables(int solverNumber) const;

  /**
   * \return The number of parameters of a solver, e.g. material values etc.
   */
  int getParameters(int solverNumber) const;

  /**
   * \return The order of the ansatz polynomials of a solver.
   */
  int getOrder(int solverNumber) const;

  /**
   * \return The maximum extent in each coordinate direction a cell is allowed
   * to have.
   */
  double getMaximumMeshSize(int solverNumber) const;

  /**
   * Prints a summary of the parameters read in for a solver.
   */
  void logSolverDetails(int solverNumber) const;

  /**
   * Checks for inconsistencies between the ExaHyPE specification file
   * and the build. Stops the program with an error
   * if both are inconsistent.
   *
   * The fields type, identifier, variables, parameters, and order
   * are considered in the inconsistency check.
   */
  void checkSolverConsistency(int solverNumber) const;

  /**
   * \return The time stepping mode of a solver.
   */
  exahype::solvers::Solver::TimeStepping getTimeStepping(
      int solverNumber) const;

  /**
   * \return The relaxation parameter used for the discrete maximum principle (DMP).
   *
   * \note This value can only be read in if the solver \p solverNumber is
   * a limiting ADER-DG solver.
   */
  double getDMPRelaxationParameter(int solverNumber) const;

  /**
   * \return The maximum-minimum difference scaling used for the discrete maximum principle (DMP).
   *
   * \note This value can only be read in if the solver \p solverNumber is
   * a limiting ADER-DG solver.
   */
  double getDMPDifferenceScaling(int solverNumber) const;

  /**
   * In the ExaHyPE specification file, a plotter configuration has
   * the following signature:
   *
   * plot <identifier> <name>
   *  variables = <variables>
   *  time      = <first-snapshot-time>
   *  repeat    = <repeat-time>
   *  output    = <filename>
   *  select    = <selector>
   * end plot
   */
  std::string getIdentifierForPlotter(int solverNumber,
                                      int plotterNumber) const;
  std::string getNameForPlotter(int solverNumber,
                                int plotterNumber) const;
  int getUnknownsForPlotter(int solverNumber, int plotterNumber) const;
  double getFirstSnapshotTimeForPlotter(int solverNumber,
                                        int plotterNumber) const;
  double getRepeatTimeForPlotter(int solverNumber, int plotterNumber) const;
  std::string getFilenameForPlotter(int solverNumber, int plotterNumber) const;
  std::string getSelectorForPlotter(int solverNumber, int plotterNumber) const;

  std::string getProfilerIdentifier() const;
  std::string getMetricsIdentifierList() const;
  std::string getProfilingOutputFilename() const;

  ParserView getParserView(int solverNumber);

  /**
   * Returns an empty string if no log file is specified in the file.
   */
  std::string getLogFileName() const;
};

#endif
