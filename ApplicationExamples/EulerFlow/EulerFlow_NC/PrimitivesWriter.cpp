#include "PrimitivesWriter.h"
#include "Primitives.h"


Euler::PrimitivesWriter::PrimitivesWriter(MyEulerSolver&  solver) {
  // @todo Please insert your code here
}


Euler::PrimitivesWriter::~PrimitivesWriter() {
  // @todo Please insert your code here
}


void Euler::PrimitivesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void Euler::PrimitivesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void Euler::PrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {

/**
 * This is the Primitive Variable plotter
 **/
  //double V[Euler::MyEulerSolver::NumberOfVariables];
  cons2prim(outputQuantities, Q);

  //for (int i=0; i<0; i++){ 
  //  outputQuantities[i] = V[i];
  //}
}


