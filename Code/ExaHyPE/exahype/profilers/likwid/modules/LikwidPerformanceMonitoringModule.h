#ifndef _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_PERFORMANCE_MONITORING_MODULE_H_
#define _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_PERFORMANCE_MONITORING_MODULE_H_

#ifdef LIKWID_AVAILABLE

#include <iostream>
#include <string>
#include <unordered_map>

#include "LikwidModule.h"

namespace exahype {
namespace profilers {
namespace likwid {

struct LikwidProfilerState;

class LikwidPerformanceMonitoringModule : public LikwidModule {
 public:
  explicit LikwidPerformanceMonitoringModule(const LikwidProfilerState& state,
                                             const std::string& group_name);
  virtual ~LikwidPerformanceMonitoringModule();

  void setNumberOfTags(int n) override;
  void registerTag(const std::string& tag) override;
  void start(const std::string& tag) override;
  void stop(const std::string& tag) override;
  void writeToOstream(std::ostream* os) const override;

 private:
  const std::string group_name_;
  int group_index_ = -1;  // index of group_name_ in groups array

  std::unordered_map<std::string, int>
      group_handles_;  // stores handles returned by addEventSet
};

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // LIKWID_AVAILABLE

#endif  // _EXAHYPE_PROFILERS_LIKWID_MODULES_LIKWID_PERFORMANCE_MONITORING_MODULE_H_