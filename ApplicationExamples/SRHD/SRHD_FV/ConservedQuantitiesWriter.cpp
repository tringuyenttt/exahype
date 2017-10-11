#include "ConservedQuantitiesWriter.h"


SRHD::ConservedQuantitiesWriter::ConservedQuantitiesWriter(SRHDSolverFV&  solver) {
  // @todo Please insert your code here
}


SRHD::ConservedQuantitiesWriter::~ConservedQuantitiesWriter() {
  // @todo Please insert your code here
}


void SRHD::ConservedQuantitiesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::ConservedQuantitiesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void SRHD::ConservedQuantitiesWriter::mapQuantities(
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


