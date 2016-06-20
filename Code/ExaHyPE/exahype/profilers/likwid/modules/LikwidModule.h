#ifndef _EXAHYPE_PROFILERS_PROFILER_LIKWID_MODULES_LIKWID_MODULE_H_
#define _EXAHYPE_PROFILERS_PROFILER_LIKWID_MODULES_LIKWID_MODULE_H_

#include <iostream>
#include <string>

namespace exahype {
namespace profilers {
namespace likwid {

struct LikwidProfilerState;

class LikwidModule {
 public:
  explicit LikwidModule(const LikwidProfilerState& state) : state_(state) {}
  virtual ~LikwidModule() {}

  // Disallow copy and assign
  LikwidModule(const LikwidModule& other) = delete;
  const LikwidModule& operator=(const LikwidModule& other) = delete;

  virtual void setNumberOfTags(int n) = 0;
  virtual void registerTag(const std::string& tag) = 0;
  virtual void start(const std::string& tag) = 0;
  virtual void stop(const std::string& tag) = 0;
  virtual void writeToOstream(std::ostream* os) const = 0;

 protected:
  const LikwidProfilerState& state_;
};

}  // namespace likwid
}  // namespace profilers
}  // namespace exahype

#endif  // _EXAHYPE_PROFILERS_PROFILER_LIKWID_MODULES_LIKWID_MODULE_H_
