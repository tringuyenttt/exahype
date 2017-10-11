#include "VariableWriter.h"
#include "AbstractCellularAutomata_FV.h"

CA::VariableWriter::VariableWriter(CellularAutomata_FV&  solver) {
  // @todo Please insert your code here
}


CA::VariableWriter::~VariableWriter() {
  // @todo Please insert your code here
}


void CA::VariableWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void CA::VariableWriter::finishPlotting() {
  // @todo Please insert your code here
}


void CA::VariableWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<(CA::AbstractCellularAutomata_FV::NumberOfVariables + CA::AbstractCellularAutomata_FV::NumberOfParameters); i++){ 
    outputQuantities[i] = Q[i];
  }
}


