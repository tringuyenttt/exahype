#include "kernels/aderdg/generic/Kernels.h"

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

extern "C" {
void adersurfaceintegral_(double* lduh, double* lFhi, double* dx);
}

void kernels::aderdg::generic::fortran::surfaceIntegral(
    double* lduh, const double* const lFbnd,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const int numberOfVariables, const int basisSize) {
    
  // circumvent 'const double'
  double* lFbndFortran = new double[numberOfVariables * 6 * basisSize * basisSize];
  for (int i = 0; i < numberOfVariables * 6 * basisSize * basisSize; i++) {
    lFbndFortran[i] = lFbnd[i];
  }

  double* dxTemp = new double[3];
  dxTemp[0] = dx[0];
  dxTemp[1] = dx[1];
  dxTemp[2] = dx[2];

  adersurfaceintegral_(lduh, lFbndFortran, dxTemp);

  delete[] lFbndFortran;
  delete[] dxTemp;
}
