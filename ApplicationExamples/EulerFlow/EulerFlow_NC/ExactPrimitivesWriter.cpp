#include "ExactPrimitivesWriter.h"
#include "Primitives.h"
#include "InitialData.h"


Euler::ExactPrimitivesWriter::ExactPrimitivesWriter(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::ExactPrimitivesWriter::~ExactPrimitivesWriter() {
  // @todo Please insert your code here
}


void Euler::ExactPrimitivesWriter::startPlotting(double time) {
  this->time = time;
}


void Euler::ExactPrimitivesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::ExactPrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  /**
   * This is the plotter for the exact solutions, given
   * as primitive Variables
   **/

  const double *xpos = x.data();
  idfunc(xpos, outputQuantities, this->time);
  cons2prim(outputQuantities, Q);

  /*
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
  */
}


