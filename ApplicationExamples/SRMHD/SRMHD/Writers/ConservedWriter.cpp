#include "Writers/ConservedWriter.h"


SRMHD::ConservedWriter::ConservedWriter(MHDSolver&  solver) {}
SRMHD::ConservedWriter::~ConservedWriter() {}
void SRMHD::ConservedWriter::startPlotting(double time) {}
void SRMHD::ConservedWriter::finishPlotting() {}

void SRMHD::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<9; i++){ 
    // Conserved quantities, nothing happens here.
    outputQuantities[i] = Q[i];
  }
}


