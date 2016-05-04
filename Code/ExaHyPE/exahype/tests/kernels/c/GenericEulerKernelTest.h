// @todo ExaHyPE Lizenz
#ifndef _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_

#include "tarch/tests/TestCase.h"
#include "peano/utils/Globals.h"

namespace exahype {
namespace tests {
namespace c {

class GenericEulerKernelTest : public tarch::tests::TestCase {
 public:
  GenericEulerKernelTest();
  virtual ~GenericEulerKernelTest();

  void run() override;

 private:
  void testPDEFluxes();
  void testSpaceTimePredictor();
  void testVolumeIntegral();
  void testRiemannSolver();
  void testSurfaceIntegral();
  void testSolutionUpdate();
#if DIMENSIONS == 2
  static void testFlux(const double* const Q, double* f, double* g);
#elif DIMENSIONS == 3
  static void testFlux(const double *const Q, double *f, double *g, double *h);
#endif

  static void testEigenvalues(const double* const Q,
                              const int normalNonZeroIndex, double* lambda);

  const double eps = 1.0e-10;  // for quick adaption of the test cases (say,
                               // switch to single precision)
};

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif