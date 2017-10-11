#include "Writers/ExactPrimitivesWriter.h"
#include "InitialDataAdapter.h"
#include "C2P-MHD.h"

SRMHD::ExactPrimitivesWriter::ExactPrimitivesWriter(MHDSolver&  solver) {}
SRMHD::ExactPrimitivesWriter::~ExactPrimitivesWriter() {}
void SRMHD::ExactPrimitivesWriter::startPlotting(double time) {}
void SRMHD::ExactPrimitivesWriter::finishPlotting() {}

void SRMHD::ExactPrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  const double *xpos = x.data();
  double Conserved[SRMHD::AbstractMHDSolver::NumberOfVariables];

  // Caveat: I don't properly use initialdatabyexahypespecfile here as it doesn't pass
  // the time. However, the Alfen Wave is the only exact solution we currently have, anyway

  alfenwave_(xpos, Conserved, &timeStamp);
  int error;
  pdecons2prim_(outputQuantities, Conserved, &error);
}


