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
#include "exahype/Parser.h"

#include "tarch/Assertions.h"

#include <fstream>

#include <stdio.h>
#include <string.h>
#include <string>
#include <regex>

#include "tarch/la/ScalarOperations.h"

tarch::logging::Log exahype::Parser::_log("exahype::Parser");

bool exahype::Parser::_interpretationErrorOccured(false);

const std::string exahype::Parser::_noTokenFound("notoken");

double exahype::Parser::getValueFromPropertyString(
    const std::string& parameterString, const std::string& key) {
  std::size_t startIndex = parameterString.find(key);
  startIndex = parameterString.find(":", startIndex);
  std::size_t endIndexBracket = parameterString.find("}", startIndex + 1);
  std::size_t endIndexComma = parameterString.find(",", startIndex + 1);

  std::size_t endIndex =
      endIndexBracket < endIndexComma ? endIndexBracket : endIndexComma;

  std::string substring =
      parameterString.substr(startIndex + 1, endIndex - startIndex - 1);

  double result;
  std::istringstream ss(substring);
  ss >> result;

  if (ss) {
    return result;
  } else {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

exahype::Parser::Parser() {
  _identifier2Type.insert(
      std::pair<std::string, exahype::solvers::Solver::Type>(
          "ADER-DG", exahype::solvers::Solver::Type::ADERDG));
  _identifier2Type.insert(
      std::pair<std::string, exahype::solvers::Solver::Type>(
          "Finite-Volumes", exahype::solvers::Solver::Type::FiniteVolumes));
  _identifier2Type.insert(
        std::pair<std::string, exahype::solvers::Solver::Type>(
            "Limiting-ADER-DG", exahype::solvers::Solver::Type::LimitingADERDG));

  _identifier2TimeStepping.insert(
      std::pair<std::string, exahype::solvers::Solver::TimeStepping>(
          "global", exahype::solvers::Solver::TimeStepping::Global));
  _identifier2TimeStepping.insert(
      std::pair<std::string, exahype::solvers::Solver::TimeStepping>(
          "globalfixed", exahype::solvers::Solver::TimeStepping::GlobalFixed));
}

void exahype::Parser::readFile(const std::string& filename) {
  const int MAX_CHARS_PER_LINE = 65536;

  std::regex COMMENT_BEGIN(R"((\/\*))"); // Covers all cases /*,/**,/***,... .
  std::regex COMMENT_END(R"((\*\/))"); //
  std::regex GROUP_BEGIN_OR_END(R"(^(\s|\t)*([a-zA-Z][^\=]+)+)");
  std::regex CONST_PARAMETER(R"(^(\s|\t)*([A-Za-z](\w|[^a-zA-Z\d\s\t])*)(\s|\t)+const(\s|\t)*=(\s|\t)*(([^\s\t]|\,\s*)+)(\s|\t)*$)");
  std::regex PARAMETER(R"(^(\s|\t)*([A-Za-z](\w|[^a-zA-Z\d\s\t])*)(\s|\t)*=(\s|\t)*(([^\s\t]|\,\s*)+)(\s|\t)*$)");
  std::regex NO_SPLITTING(R"(\}|\{)");
  std::regex COMMA_SEPARATED(R"((\w|[^a-zA-Z\,\s\t])+)");
  std::regex WHITESPACE_SEPARATED(R"(([^\s\t]+))");
  std::smatch match;

  _tokenStream.clear();

  std::ifstream inputFile;
  inputFile.open(filename.c_str());
  if (!inputFile.good()) {
    logError("readFile(String)", "cannot open file " << filename);
    _tokenStream.clear();
    _interpretationErrorOccured = true;
    return;
  }

  int currentlyReadsComment = 0;
  int lineNumber            = 0;
  while (!inputFile.eof() && inputFile) {
    char lineBuffer[MAX_CHARS_PER_LINE];
    inputFile.getline(lineBuffer, MAX_CHARS_PER_LINE);
    std::string line(lineBuffer);

    // parse the line
    if (std::regex_search(line, match, COMMENT_BEGIN) && match.size() > 1) {
      currentlyReadsComment += 1;
    }
    if (std::regex_search(line, match, COMMENT_END) && match.size() > 1) {
      currentlyReadsComment -= 1;
    }

    // Runtime parameters
    if (currentlyReadsComment==0 && std::regex_search(line, match, PARAMETER) && match.size() > 1) {
      _tokenStream.push_back(match.str(2)); // Subgroup 2 is left-hand side (trimmed)
      std::string rightHandSide = match.str(6);

      if (!std::regex_search(rightHandSide, match, NO_SPLITTING)) {
        std::regex_iterator<std::string::iterator> regex_it ( rightHandSide.begin(), rightHandSide.end(), COMMA_SEPARATED );
        std::regex_iterator<std::string::iterator> rend;
        while(regex_it!=rend) {
          _tokenStream.push_back(regex_it->str());
          ++regex_it;
        }
      } else {
        _tokenStream.push_back(rightHandSide);
      }
    // Compile time parameters (Do not push the token const on the stream)
    } else if (currentlyReadsComment==0 && std::regex_search(line, match, CONST_PARAMETER) && match.size() > 1) {
      _tokenStream.push_back(match.str(2)); // Subgroup 2 is left-hand side (trimmed)
      std::string rightHandSide = match.str(7);

      if (!std::regex_search(rightHandSide, match, NO_SPLITTING)) {
        std::regex_iterator<std::string::iterator> regex_it ( rightHandSide.begin(), rightHandSide.end(), COMMA_SEPARATED );
        std::regex_iterator<std::string::iterator> rend;
        while(regex_it!=rend) {
          _tokenStream.push_back(regex_it->str());
          ++regex_it;
        }
      } else {
        _tokenStream.push_back(rightHandSide);
      }
    } else if (currentlyReadsComment==0 && std::regex_search(line, match, GROUP_BEGIN_OR_END) && match.size() > 1) {
      std::regex_iterator<std::string::iterator> regex_it ( line.begin(), line.end(), WHITESPACE_SEPARATED );
      std::regex_iterator<std::string::iterator> rend;
      if (regex_it->str().compare("end")!=0) { // first token should not be end
        while(regex_it!=rend) {
          _tokenStream.push_back(regex_it->str());
          ++regex_it;
        }
      } // else do nothing
    } else if (currentlyReadsComment<0) {
      logError("readFile(String)",
           "Please remove additional multi-line comment end(s) in line '" << lineNumber << "'.");
       _interpretationErrorOccured = true;
    }
    lineNumber++;
  }

  if (currentlyReadsComment>0) {
    logError("readFile(String)",
             "A multi-line comment was not closed after line " << lineNumber);
    _interpretationErrorOccured = true;
  }

//  For debugging purposes
/*
  std::cout << "_tokenStream=" << std::endl;
  for (std::string str : _tokenStream) {
    std::cout << "["<<str<<"]" << std::endl;
  }
*/
}

bool exahype::Parser::isValid() const {
  return !_tokenStream.empty() && !_interpretationErrorOccured;
}

std::string exahype::Parser::getTokenAfter(std::string token,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::Parser::getTokenAfter(std::string token0,
                                           std::string token1,
                                           int additionalTokensToSkip) const {
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token0) {
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         _tokenStream[currentToken] != token1) {
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::Parser::getTokenAfter(std::string token0, int occurance0,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

std::string exahype::Parser::getTokenAfter(std::string token0, int occurance0,
                                           std::string token1, int occurance1,
                                           int additionalTokensToSkip) const {
  assertion(isValid());
  assertion(occurance0 > 0);
  assertion(occurance1 > 0);
  int currentToken = 0;
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token0 || occurance0 > 1)) {
    if (_tokenStream[currentToken] == token0) occurance0--;
    currentToken++;
  }
  while (currentToken < static_cast<int>(_tokenStream.size()) &&
         (_tokenStream[currentToken] != token1 || occurance1 > 1)) {
    if (_tokenStream[currentToken] == token1) occurance1--;
    currentToken++;
  }
  currentToken += (additionalTokensToSkip + 1);
  if (currentToken < static_cast<int>(_tokenStream.size())) {
    return _tokenStream[currentToken];
  } else
    return _noTokenFound;
}

int exahype::Parser::getNumberOfThreads() const {
  assertion(isValid());
  std::string token = getTokenAfter("shared-memory", "cores");
  logDebug("getNumberOfThreads()", "found token " << token);
  int result = atoi(token.c_str());
  if (result == 0) {
    logError("getNumberOfThreads()",
             "Invalid number of cores set or token shared-memory missing: "
                 << token);
    result = 1;
    _interpretationErrorOccured = true;
  }
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getDomainSize() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;
  token = getTokenAfter("computational-domain", "width", 0);
  result(0) = atof(token.c_str());
  token = getTokenAfter("computational-domain", "width", 1);
  result(1) = atof(token.c_str());
#if DIMENSIONS == 3
  token = getTokenAfter("computational-domain", "width", 2);
  result(2) = atof(token.c_str());
#endif
  return result;
}

tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getBoundingBoxSize()
    const {
  tarch::la::Vector<DIMENSIONS, double> domainSize = getDomainSize();
  double longestH = tarch::la::max(domainSize);
  return tarch::la::Vector<DIMENSIONS, double>(longestH);
}

tarch::la::Vector<DIMENSIONS, double> exahype::Parser::getOffset() const {
  assertion(isValid());
  std::string token;
  tarch::la::Vector<DIMENSIONS, double> result;
  token = getTokenAfter("computational-domain", "offset", 0);
  result(0) = atof(token.c_str());
  token = getTokenAfter("computational-domain", "offset", 1);
  result(1) = atof(token.c_str());
#if DIMENSIONS == 3
  token = getTokenAfter("computational-domain", "offset", 2);
  result(2) = atof(token.c_str());
#endif
  logDebug("getSize()", "found offset " << result);
  return result;
}

std::string exahype::Parser::getMulticorePropertiesFile() const {
  std::string result = getTokenAfter("shared-memory", "properties-file");
  logDebug("getMulticorePropertiesFile()", "found token " << result);
  return result;
}

exahype::Parser::MPILoadBalancingType exahype::Parser::getMPILoadBalancingType()
    const {
  std::string token = getTokenAfter("distributed-memory", "identifier");
  exahype::Parser::MPILoadBalancingType result = MPILoadBalancingType::Static;
  if (token.compare("static_load_balancing") == 0) {
    result = MPILoadBalancingType::Static;
  } else {
    logError("getMPILoadBalancingType()",
             "Invalid distributed memory identifier " << token);
    _interpretationErrorOccured = true;
  }
  return result;
}

std::string exahype::Parser::getMPIConfiguration() const {
  return getTokenAfter("distributed-memory", "configure");
}

int exahype::Parser::getMPIBufferSize() const {
  std::string token = getTokenAfter("distributed-memory", "buffer-size");
  int result = atoi(token.c_str());
  if (result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI buffer size " << token);
    result = 64;
    _interpretationErrorOccured = true;
  }
  return result;
}

int exahype::Parser::getMPITimeOut() const {
  std::string token = getTokenAfter("distributed-memory", "timeout");
  int result = atoi(token.c_str());
  if (result <= 0) {
    logError("getMPIBufferSize()", "Invalid MPI timeout value " << token);
    result = 0;
    _interpretationErrorOccured = true;
  }
  return result;
}

exahype::Parser::MulticoreOracleType exahype::Parser::getMulticoreOracleType()
    const {
  std::string token = getTokenAfter("shared-memory", "identifier");
  exahype::Parser::MulticoreOracleType result = MulticoreOracleType::Dummy;
  if (token.compare("dummy") == 0) {
    result = MulticoreOracleType::Dummy;
  } else if (token.compare("autotuning") == 0) {
    result = MulticoreOracleType::Autotuning;
  } else if (token.compare("sampling") == 0) {
    result = MulticoreOracleType::GrainSizeSampling;
  } else {
    logError("getMulticoreOracleType()", "Invalid shared memory identifier "
                                             << token);
    result = MulticoreOracleType::Dummy;
    _interpretationErrorOccured = true;
  }
  return result;
}

double exahype::Parser::getSimulationEndTime() const {
  std::string token = getTokenAfter("computational-domain", "end-time");
  logDebug("getSimulationEndTime()", "found token " << token);
  double result = atof(token.c_str());
  if (result <= 0) {
    logError("getSimulationEndTime()",
             "Invalid simulation end-time: " << token);
    result = 1.0;
    _interpretationErrorOccured = true;
  }
  return result;
}

bool exahype::Parser::getFuseAlgorithmicSteps() const {
  std::string token = getTokenAfter("optimisation", "fuse-algorithmic-steps");
  logDebug("getFuseAlgorithmicSteps()", "found fuse-algorithmic-steps"
                                            << token);

  bool result = token.compare("on") == 0;

  if (token.compare(_noTokenFound) == 0) {
    result = false;  // default value
  } else if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getFuseAlgorithmicSteps()",
             "fuse-algorithmic-steps is required in the "
             "optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }
  return result;
}

bool exahype::Parser::getExchangeBoundaryDataInBatchedTimeSteps() const {
  std::string token = getTokenAfter(
      "optimisation",
      "disable-amr-if-grid-has-been-stationary-in-previous-iteration");
  if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getExchangeBoundaryDataInBatchedTimeSteps()",
             "disable-amr-if-grid-has-been-stationary-in-previous-iteration is "
             "required in the "
             "optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }
  return token.compare("off") == 0;
}

double exahype::Parser::getFuseAlgorithmicStepsFactor() const {
  std::string token =
      getTokenAfter("optimisation", "fuse-algorithmic-steps-factor");

  char* pEnd;
  double result = std::strtod(token.c_str(), &pEnd);
  logDebug("getFuseAlgorithmicStepsFactor()",
           "found fuse-algorithmic-steps-factor " << token);

  if (result < 0.0 || result > 1.0 || pEnd == token.c_str()) {
    logError("getFuseAlgorithmicStepsFactor()",
             "'fuse-algorithmic-steps-factor': Value must be greater than zero "
             "and smaller than one: "
                 << result);
    result = 0.0;
    _interpretationErrorOccured = true;
  }

  return result;
}

double exahype::Parser::getTimestepBatchFactor() const {
  std::string token = getTokenAfter("optimisation", "timestep-batch-factor");
  char* pEnd;
  double result = std::strtod(token.c_str(), &pEnd);
  logDebug("getFuseAlgorithmicStepsFactor()", "found timestep-batch-factor "
                                                  << token);

  if (result < 0.0 || result > 1.0 || pEnd == token.c_str()) {
    logError("getFuseAlgorithmicStepsFactor()",
             "'timestep-batch-factor': Value is required in optimisation "
             "section and must be greater than zero and smaller than one: "
                 << result);
    result = 0.0;
    _interpretationErrorOccured = true;
  }

  return result;
}

bool exahype::Parser::getSkipReductionInBatchedTimeSteps() const {
  std::string token =
      getTokenAfter("optimisation", "skip-reduction-in-batched-time-steps");
  logDebug("getSkipReductionInBatchedTimeSteps()",
           "found skip-reduction-in-batched-time-steps " << token);
  if (token.compare("on") != 0 && token.compare("off") != 0) {
    logError("getSkipReductionInBatchedTimeSteps()",
             "skip-reduction-in-batched-time-steps is required in the "
             "optimisation segment and has to be either on or off: "
                 << token);
    _interpretationErrorOccured = true;
  }

  return token.compare("on") == 0;
}


double exahype::Parser::getDoubleCompressionFactor() const {
  std::string token = getTokenAfter("optimisation", "double-compression");

  if (token.compare(_noTokenFound) == 0) {
    return 0.0;  // default value
  }
  else {
    char* pEnd;
    double result = std::strtod(token.c_str(), &pEnd);
    logDebug("getDoubleCompressionFactor()", "found double-compression "
                                                  << token);

    if (result < 0.0 || pEnd == token.c_str()) {
      logError("getDoubleCompressionFactor()",
             "'double-compression': Value is required in optimisation "
             "section and must be greater than or equal to zero: " << result);
      result = 0.0;
      _interpretationErrorOccured = true;
    }

    return result;
  }
}


bool   exahype::Parser::getSpawnDoubleCompressionAsBackgroundTask() const {
  std::string token =
      getTokenAfter("optimisation", "spawn-double-compression-as-background-thread");

  if (token.compare(_noTokenFound) == 0) {
    return false;  // default value
  }
  else {
    logDebug("getSpawnDoubleCompressionAsBackgroundTask()",
           "found spawn-double-compression-as-background-thread " << token);
    if (token.compare("on") != 0 && token.compare("off") != 0) {
      logError("getSpawnDoubleCompressionAsBackgroundTask()",
             "spawn-double-compression-as-background-thread is required in the "
             "optimisation segment and has to be either on or off: "
                 << token);
      _interpretationErrorOccured = true;
    }

    return token.compare("on") == 0;
  }
}


exahype::solvers::Solver::Type exahype::Parser::getType(
    int solverNumber) const {
  std::string token;
  exahype::solvers::Solver::Type result =
      exahype::solvers::Solver::Type::ADERDG;
  token = getTokenAfter("solver", solverNumber + 1, 0);
  if (_identifier2Type.find(token) != _identifier2Type.end()) {
    result = _identifier2Type.at(token);
    // logDebug("getType()", "found type " << result);
    logDebug("getType()", "found type ");
  } else {
    logError(
        "getType()",
        "'" << getIdentifier(solverNumber) << "': 'type': Value '" << token
            << "' is invalid. See the ExaHyPE documentation for valid values.");
    _interpretationErrorOccured = true;
  }
  return result;
}

std::string exahype::Parser::getIdentifier(int solverNumber) const {
  std::string token;
  token = getTokenAfter("solver", solverNumber + 1, 1);
  logDebug("getIdentifier()", "found identifier " << token);
  return token;
}

int exahype::Parser::getVariables(int solverNumber) const {
  std::string token;
  int result;
  std::regex COLON_SEPARATED(R"(([A-Za-z]\w*):([0-9]+))");
  std::smatch match;

  // first check if we read in a number
  token = getTokenAfter("solver", solverNumber + 1, "variables", 1);
  result = atoi(token.c_str());

  if (result < 1) { // token is not a number
    result = 0;

    int i = 1;
    std::regex_search(token, match, COLON_SEPARATED);
    while (match.size() > 1) {
      int multiplicity = atoi(match.str(2).c_str());
      result +=multiplicity; // std::string name = match.str(1);

      // logInfo("getVariables(...)","token="<<token<<",name="<<match.str(1)<<",n="<<match.str(2));
      token = getTokenAfter("solver", solverNumber + 1, "variables", 1, i++);
      std::regex_search(token, match, COLON_SEPARATED);
    }

    if (result < 1) { // token is still 0
      logError("getVariables()",
               "'" << getIdentifier(solverNumber)
               << "': 'variables': Value must be greater than zero.");
          _interpretationErrorOccured = true;
    }
  }

  logDebug("getVariables()", "found variables " << result);
  return result;
}

int exahype::Parser::getParameters(int solverNumber) const {
  std::string token;
  int result;
  std::regex COLON_SEPARATED(R"(([A-Za-z]\w*):([0-9]+))");
  std::smatch match;

  // first check if we read in a number
  token = getTokenAfter("solver", solverNumber + 1, "parameters", 1);
  result = atoi(token.c_str());

  if (result < 1) { // token is not a number
    result = 0;

    int i = 1;
    std::regex_search(token, match, COLON_SEPARATED);
    while (match.size() > 1) {
      int multiplicity = atoi(match.str(2).c_str());
      result +=multiplicity; // std::string name = match.str(1);

      // logInfo("getVariables(...)","token="<<token<<",name="<<match.str(1)<<",n="<<match.str(2));
      token = getTokenAfter("solver", solverNumber + 1, "parameters", 1, i++);
      std::regex_search(token, match, COLON_SEPARATED);
    }

    if (result < 0) { // token is still 0
      logError("getParameters()",
               "'" << getIdentifier(solverNumber)
               << "': 'parameters': Value must be non-negative.");
          _interpretationErrorOccured = true;
    }
  }

  logDebug("getVariables()", "found variables " << result);
  return result;
}

int exahype::Parser::getOrder(int solverNumber) const {
  std::string token;
  int result;
  token = getTokenAfter("solver", solverNumber + 1, "order", 1);
  result = atoi(token.c_str());

  if (result < 0) {
    logError("getOrder()", "'" << getIdentifier(solverNumber)
                               << "': 'order': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logDebug("getOrder()", "found order " << result);
  return result;
}


double exahype::Parser::getMaximumMeshSize(int solverNumber) const {
  std::string token;
  double result;
  token =
      getTokenAfter("solver", solverNumber + 1, "maximum-mesh-size", 1, 0);
  result = atof(token.c_str());
  if (tarch::la::smallerEquals(result, 0.0)) {
    logError("getMaximumMeshSize()",
             "'" << getIdentifier(solverNumber)
                 << "': 'maximum-mesh-size': Value must be greater than zero.");
    _interpretationErrorOccured = true;
  }

  logDebug("getMaximumMeshSize()", "found maximum mesh size " << result);
  return result;
}

exahype::solvers::Solver::TimeStepping exahype::Parser::getTimeStepping(
    int solverNumber) const {
  std::string token;
  exahype::solvers::Solver::TimeStepping result;
  token = getTokenAfter("solver", solverNumber + 1, "time-stepping", 1);
  if (_identifier2TimeStepping.find(token) != _identifier2TimeStepping.end()) {
    result = _identifier2TimeStepping.at(token);
    // logDebug("getTimeStepping()", "found TimeStepping " << result);
    logDebug("getTimeStepping()", "found TimeStepping "<< token);
    return result;
  } else {
    logError(
        "getTimeStepping()",
        "'" << getIdentifier(solverNumber) << "': 'time-stepping': Value '"
            << token
            << "' is invalid. See the ExaHyPE documentation for valid values.");
    _interpretationErrorOccured = true;
  }
  return exahype::solvers::Solver::TimeStepping::Global;
}

double exahype::Parser::getDMPRelaxationParameter(int solverNumber) const {
  std::string token;
  double result;
  token = getTokenAfter("solver", solverNumber + 1, "dmp-relaxation-parameter", 1);
  result = atof(token.c_str());

  if (result < 0) {
    logError("getDMPRelaxationParameter()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-relaxation-parameter': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logInfo("getParameters()", "found dmp-relaxation-parameter " << result);
  return result;
}

double exahype::Parser::getDMPDifferenceScaling(int solverNumber) const {
  std::string token;
  double result;
  token = getTokenAfter("solver", solverNumber + 1, "dmp-difference-scaling", 1);
  result = atof(token.c_str());

  if (result < 0) {
    logError("getParameters()",
             "'" << getIdentifier(solverNumber)
                 << "': 'dmp-difference-scaling': Value must not be negative.");
    _interpretationErrorOccured = true;
  }

  logInfo("getDMPDifferenceScaling()", "found dmp-difference-scaling " << result);
  return result;
}

std::string exahype::Parser::getIdentifierForPlotter(int solverNumber,
                                                     int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1);
  logDebug("getIdentifierForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

std::string exahype::Parser::getNameForPlotter(int solverNumber,
                                               int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 1);
  logDebug("getIdentifierForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

int exahype::Parser::getUnknownsForPlotter(int solverNumber,
                                           int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 3);
  logDebug("getUnknownsForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return atoi(token.c_str());
}

double exahype::Parser::getFirstSnapshotTimeForPlotter(
    int solverNumber, int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 5);
  logDebug("getFirstSnapshotTimeForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return atof(token.c_str());
}

double exahype::Parser::getRepeatTimeForPlotter(int solverNumber,
                                                int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 7);
  logDebug("getRepeatTimeForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return atof(token.c_str());
}

std::string exahype::Parser::getFilenameForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 9);
  logDebug("getFilenameForPlotter()", "found token " << token);
  assertion3(token.compare(_noTokenFound) != 0, token, solverNumber,
             plotterNumber);
  return token;
}

std::string exahype::Parser::getSelectorForPlotter(int solverNumber,
                                                   int plotterNumber) const {
  // We have to multiply with two as the token solver occurs twice (to open and
  // close the section)
  std::string token = getTokenAfter("solver", solverNumber + 1, "plot",
                                    plotterNumber + 1, 11);
  logDebug("getSelectorForPlotter()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

std::string exahype::Parser::getLogFileName() const {
  std::string token = getTokenAfter("log-file");
  logDebug("getLogFileName()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

std::string exahype::Parser::getProfilerIdentifier() const {
  std::string token = getTokenAfter("profiling", "profiler");
  logDebug("getProfilerIdentifier()", "found token " << token);
  return (token != _noTokenFound) ? token : "NoOpProfiler";
}

std::string exahype::Parser::getMetricsIdentifierList() const {
  std::string token = getTokenAfter("profiling", "metrics");
  logDebug("getMetricsIdentifierList()", "found token " << token);
  return (token != _noTokenFound) ? token : "{}";
}

std::string exahype::Parser::getProfilingOutputFilename() const {
  std::string token = getTokenAfter("profiling", "profiling-output");
  logDebug("getProfilingOutputFilename()", "found token " << token);
  return (token != _noTokenFound) ? token : "";
}

void exahype::Parser::logSolverDetails(int solverNumber) const {
  logInfo("logSolverDetails()",
          "Solver " << getTokenAfter("solver", solverNumber + 1, 0) << " "
                    << getIdentifier(solverNumber) << ":");
  logInfo("logSolverDetails()", "variables:\t\t" << getVariables(solverNumber));
  logInfo("logSolverDetails()", "parameters:\t" << getParameters(solverNumber));
  logInfo("logSolverDetails()", "order:\t\t" << getOrder(solverNumber));
  logInfo("logSolverDetails()", "maximum-mesh-size:\t"
                                    << getMaximumMeshSize(solverNumber));
  logInfo("logSolverDetails()",
          "time-stepping:\t" << getTokenAfter("solver", solverNumber + 1,
                                              "time-stepping", 1));
}

void exahype::Parser::checkSolverConsistency(int solverNumber) const {
  assertion1(solverNumber <
                 static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
             solverNumber);
  exahype::solvers::Solver* solver =
      exahype::solvers::RegisteredSolvers[solverNumber];

  bool recompile = false;
  bool runToolkitAgain = false;
  if (solver->getType() != getType(solverNumber)) {
    logError("checkIfSolverIsConsistent",
             "'" << getIdentifier(solverNumber)
                 << "': Solver type in specification file"
                 << "differs from implementation solver type.");
    recompile = true;
    _interpretationErrorOccured = true;
  }

  if (solver->getIdentifier().compare(getIdentifier(solverNumber))) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Identifier in specification file "
                 << "('" << getIdentifier(solverNumber)
                 << "') differs from identifier used in implementation ('"
                 << solver->getIdentifier() << "').");
    recompile = true;
    _interpretationErrorOccured = true;
  }

  if (solver->getNumberOfVariables() != getVariables(solverNumber)) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for 'variables' in specification file"
                 << "('" << getVariables(solverNumber)
                 << "') differs from number of variables used in "
                    "implementation file ('"
                 << solver->getNumberOfVariables() << "').");
    recompile = true;
    _interpretationErrorOccured = true;
  }

  if (solver->getNumberOfParameters() != getParameters(solverNumber)) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for field 'parameters' in specification file"
                 << "('" << getParameters(solverNumber)
                 << "') differs from  number of parameters used in "
                    "implementation file ('"
                 << solver->getNumberOfParameters() << "').");
    recompile = true;
    _interpretationErrorOccured = true;
  }

  if (solver->getType() == exahype::solvers::Solver::Type::ADERDG &&
      solver->getNodesPerCoordinateAxis() != getOrder(solverNumber) + 1) {
    logError("checkSolverConsistency",
             "'" << getIdentifier(solverNumber)
                 << "': Value for field 'order' in specification file "
                 << "('" << getOrder(solverNumber)
                 << "') differs from value used in implementation file ('"
                 << solver->getNodesPerCoordinateAxis() - 1 << "'). ");
    runToolkitAgain = true;
    _interpretationErrorOccured = true;
  }

  // @todo We should add checks for FV as well

  if (runToolkitAgain) {
    logError("checkSolverConsistency",
             "Please (1) run the Toolkit again, and (2) recompile!");
    _interpretationErrorOccured = true;
  }

  if (recompile) {
    logError(
        "checkSolverConsistency",
        "Please (1) adjust the specification file (*.exahype) or the file '"
            << solver->getIdentifier()
            << ".cpp' accordingly, and (2) recompile!");
    _interpretationErrorOccured = true;
  }
}

exahype::Parser::ParserView exahype::Parser::getParserView(int solverNumber) {
  return ParserView(*this, solverNumber);
}

exahype::Parser::ParserView::ParserView(Parser& parser,
                                        int solverNumberInSpecificationFile)
    : _parser(parser),
      _solverNumberInSpecificationFile(solverNumberInSpecificationFile) {}

std::string exahype::Parser::ParserView::getValue(
    const std::string inputString, const std::string& key) const {
  assertion(_parser.isValid());

  if (inputString.substr(0, 1) != "{") return "";
  std::size_t currentIndex = 1;
  bool nextTokenIsSearchedValue = false;
  while (currentIndex < inputString.size() - 1) {
    std::size_t firstIndexColon = inputString.find(":", currentIndex);
    std::size_t firstIndexBracket = inputString.find("}", currentIndex);
    std::size_t firstIndexComma = inputString.find(",", currentIndex);
    std::size_t endIndex =
        std::min(firstIndexColon, std::min(firstIndexBracket, firstIndexComma));

    std::string newToken =
        inputString.substr(currentIndex, endIndex - currentIndex);
    if (nextTokenIsSearchedValue) {
      return newToken;
    }

    nextTokenIsSearchedValue =
        (endIndex == firstIndexColon) && (newToken == key);

    currentIndex = endIndex + 1;
  }

  return "";
}

bool exahype::Parser::ParserView::hasKey(const std::string& key) const {
  assertion(_parser.isValid());

  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);

  if (inputString.substr(0, 1) != "{") return false;
  std::size_t currentIndex = 1;
  bool nextTokenIsValue = false;
  while (currentIndex < inputString.size() - 1) {
    std::size_t firstIndexColon = inputString.find(":", currentIndex);
    std::size_t firstIndexBracket = inputString.find("}", currentIndex);
    std::size_t firstIndexComma = inputString.find(",", currentIndex);
    std::size_t endIndex =
        std::min(firstIndexColon, std::min(firstIndexBracket, firstIndexComma));

    if (!nextTokenIsValue) {
      std::string newToken =
          inputString.substr(currentIndex, endIndex - currentIndex);
      logDebug("hasKey(string)", "added token \"" << newToken << "\"");
      if (newToken == key) {
        return true;
      }
    }

    nextTokenIsValue = (endIndex == firstIndexColon);

    currentIndex = endIndex + 1;
  }

  return false;
}

int exahype::Parser::ParserView::getValueAsInt(const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  int result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    assertion(!isValueValidInt(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return -1;
  }
}

bool exahype::Parser::ParserView::getValueAsBool(const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  bool result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    assertion(!isValueValidBool(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return -1;
  }
}

double exahype::Parser::ParserView::getValueAsDouble(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  double result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return result;
  } else {
    assertion(!isValueValidDouble(key));
    assertionMsg(false, "shall not happen. Please call isValueValidXXX before");
    return -1.0;
  }
}

std::string exahype::Parser::ParserView::getValueAsString(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  return getValue(inputString, key);
}

bool exahype::Parser::ParserView::isValueValidInt(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  int result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    return false;
  }
}

bool exahype::Parser::ParserView::isValueValidDouble(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  double result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    return false;
  }
}


bool exahype::Parser::ParserView::isValueValidBool(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  std::string value = getValue(inputString, key);

  bool result;
  std::istringstream ss(value);
  ss >> result;

  if (ss) {
    return true;
  } else {
    return false;
  }
}

bool exahype::Parser::ParserView::isValueValidString(
    const std::string& key) const {
  const std::string inputString = _parser.getTokenAfter(
      "solver", _solverNumberInSpecificationFile + 1, "constants", 1);
  return getValue(inputString, key) != "";
}
