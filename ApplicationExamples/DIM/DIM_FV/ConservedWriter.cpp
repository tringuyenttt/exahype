#include "ConservedWriter.h"


DIM::ConservedWriter::ConservedWriter(DIMSolver&  solver) {
  // @todo Please insert your code here
}


DIM::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void DIM::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void DIM::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void DIM::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<14; i++){ 
    outputQuantities[i] = Q[i];
  }
}


