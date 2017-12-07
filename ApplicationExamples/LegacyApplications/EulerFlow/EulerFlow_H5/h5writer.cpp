#include "h5writer.h"


Euler::h5writer::h5writer(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::h5writer::~h5writer() {
  // @todo Please insert your code here
}


void Euler::h5writer::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::h5writer::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::h5writer::mapQuantities(
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
  
  // DEBUGGING:
  // abuse some output quantities:
  const int start=4;
  outputQuantities[start+1] = x(0);
  outputQuantities[start+2] = x(1);
  outputQuantities[start+3] = DIMENSIONS==3 ? x(2) : -1;
  outputQuantities[start+4] = pos(0);
  outputQuantities[start+5] = pos(1);
  outputQuantities[start+6] = DIMENSIONS==3 ? pos(2) : -1;
  outputQuantities[start+7] = offsetOfPatch(0);
  outputQuantities[start+8] = offsetOfPatch(1);
  outputQuantities[start+9] = DIMENSIONS==3 ? offsetOfPatch(2) : -1;
}


