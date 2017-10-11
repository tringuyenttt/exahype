#include "ConservedWriter.h"


SWE::ConservedWriter::ConservedWriter(MySWESolver&  solver) {
  // @todo Please insert your code here
}


SWE::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void SWE::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<4; i++){
    outputQuantities[i] = Q[i];
  }
  outputQuantities[4] = Q[3] + Q[0];
}


