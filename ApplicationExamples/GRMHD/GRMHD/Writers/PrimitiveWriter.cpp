#include "PrimitiveWriter.h"
#include "Fortran/PDE.h"


GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_FV&  solver) {
  // @todo Please insert your code here
}



GRMHD::PrimitiveWriter::PrimitiveWriter(GRMHDSolver_ADERDG&  solver) {
  // @todo Please insert your code here
}


GRMHD::PrimitiveWriter::PrimitiveWriter(exahype::solvers::LimitingADERDGSolver&  solver) {
  // @todo Please insert your code here
}

GRMHD::PrimitiveWriter::~PrimitiveWriter() {
  // @todo Please insert your code here
}


void GRMHD::PrimitiveWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void GRMHD::PrimitiveWriter::finishPlotting() {
  // @todo Please insert your code here
}


void GRMHD::PrimitiveWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
    int err;
    pdecons2prim_(outputQuantities, Q, &err);
    if(err != 0) {
	    printf("Cons2Prim Failure in PrimitiveWriter!!!");
    }
}


