#include "FVM_Plotter0.h"


Euler::FVM_Plotter0::FVM_Plotter0(FVM&  solver) {
  // @todo Please insert your code here
}


Euler::FVM_Plotter0::~FVM_Plotter0() {
  // @todo Please insert your code here
}


void Euler::FVM_Plotter0::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::FVM_Plotter0::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::FVM_Plotter0::mapQuantities(
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


