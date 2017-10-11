#include "ErrorWriter.h"

#include "MHDSolver.h"
#include "InitialDataAdapter.h"


SRMHD::ErrorWriter::ErrorWriter(MHDSolver&  solver) {
}


SRMHD::ErrorWriter::~ErrorWriter() {
}


void SRMHD::ErrorWriter::startPlotting(double time) {
}


void SRMHD::ErrorWriter::finishPlotting() {
}


void SRMHD::ErrorWriter::mapQuantities(
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
  alfenwave_(xpos, outputQuantities, &timeStamp);
}


