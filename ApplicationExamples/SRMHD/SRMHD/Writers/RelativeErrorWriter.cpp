#include "Writers/RelativeErrorWriter.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"

#include "C2P-MHD.h"

#include <limits> // double max
#include <cmath> // C++11: std::signbit

SRMHD::RelativeErrorWriter::RelativeErrorWriter(MHDSolver&  solver) {}
SRMHD::RelativeErrorWriter::~RelativeErrorWriter() {}
void SRMHD::RelativeErrorWriter::startPlotting(double time) {}
void SRMHD::RelativeErrorWriter::finishPlotting() {}

// There is also http://stackoverflow.com/a/25909704 not relying on C++11
template<typename T>
T sign(T arg) {
  // Returns the signum of arg. Works especially with NaNs and zeros.
  // signbit(arg) is true if arg is negative
  return std::signbit(arg) ? T(-1) : T(+1);
}

void SRMHD::RelativeErrorWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const int nDim = DIMENSIONS;
  const int nVar = SRMHD::AbstractMHDSolver::NumberOfVariables;
  
  const double *xpos = x.data();
  double exactCons[nVar];

  alfenwave_(xpos, exactCons, &timeStamp);

  // compute the difference
  for (int i=0; i<nVar; i++) {
    outputQuantities[i] = ( Q[i] - exactCons[i] ) / exactCons[i];
    // Make NaN to finite value. Basically works like numpy.nan_to_num.
    if(outputQuantities[i] != outputQuantities[i] || outputQuantities[i] == std::numeric_limits<double>::infinity()) {
        outputQuantities[i] = sign(outputQuantities[i]) * std::numeric_limits<double>::max();
    }
  }
}


