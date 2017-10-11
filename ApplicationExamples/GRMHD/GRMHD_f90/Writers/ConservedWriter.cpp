#include "ConservedWriter.h"


GRMHD::ConservedWriter::ConservedWriter(GRMHDSolver&  solver) {
  // @todo Please insert your code here
}


GRMHD::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void GRMHD::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<19; i++){ 
    outputQuantities[i] = Q[i];
  }
}


