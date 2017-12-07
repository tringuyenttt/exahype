#include "ConservedWriter.h"


Euler::ConservedWriter::ConservedWriter(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::ConservedWriter::~ConservedWriter() {
  // @todo Please insert your code here
}


void Euler::ConservedWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::ConservedWriter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
}


