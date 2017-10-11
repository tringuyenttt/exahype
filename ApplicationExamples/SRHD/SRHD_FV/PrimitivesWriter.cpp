#include "PrimitivesWriter.h"


SRHD::PrimitivesWriter::PrimitivesWriter(SRHDSolverFV&  solver) {
  // @todo Please insert your code here
}


SRHD::PrimitivesWriter::~PrimitivesWriter() {
  // @todo Please insert your code here
}


void SRHD::PrimitivesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void SRHD::PrimitivesWriter::finishPlotting() {
  // @todo Please insert your code here
}

extern "C" {
  // Fortran function, implemented in C2P-SRHD.f90                                                                                                                                    
  void pdecons2prim_(double* V, double* Q, int* iErr);
}

void SRHD::PrimitivesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  // SRHD_FV primitive variables                                                                                                           
  int iErr;
  pdecons2prim_(outputQuantities, Q, &iErr);

}


