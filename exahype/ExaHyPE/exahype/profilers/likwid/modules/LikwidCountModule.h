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

#ifndef _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_COUNT_MODULE_H_
#define _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_COUNT_MODULE_H_

#ifdef LIKWID_AVAILABLE

#include <iostream>
#include <string>
#include <unordered_map>

#include "LikwidModule.h"

namespace exahype {
namespace profilers {
namespace likwid {

class LikwidCountModule : public LikwidModule {
 public:
  explicit LikwidCountModule(const LikwidProfilerState& state)
      : LikwidModule(state) {}
  virtual ~LikwidCountModule() {}

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  std::unordered_map<std::string, unsigned long> counts_;
};

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_COUNT_MODULE_H_
