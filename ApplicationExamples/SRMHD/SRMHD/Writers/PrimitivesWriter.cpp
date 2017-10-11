#include "Writers/PrimitivesWriter.h"
#include "C2P-MHD.h"

SRMHD::PrimitivesWriter::PrimitivesWriter(MHDSolver&  solver) {}
SRMHD::PrimitivesWriter::~PrimitivesWriter() {}
void SRMHD::PrimitivesWriter::startPlotting(double time) {}
void SRMHD::PrimitivesWriter::finishPlotting() {}

void SRMHD::PrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // MHD primitive variable plotter
  int i;
  pdecons2prim_(outputQuantities, Q, &i);
}


