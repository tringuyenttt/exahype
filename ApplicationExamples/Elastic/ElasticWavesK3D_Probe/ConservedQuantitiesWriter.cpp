#include "ConservedQuantitiesWriter.h"

#include "MyElasticWaveSolver.h"

ElasticWave::ConservedQuantitiesWriter::ConservedQuantitiesWriter(MyElasticWaveSolver&  solver) {
  // @todo Please insert your code here
}


ElasticWave::ConservedQuantitiesWriter::~ConservedQuantitiesWriter() {
  // @todo Please insert your code here
}


void ElasticWave::ConservedQuantitiesWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void ElasticWave::ConservedQuantitiesWriter::finishPlotting() {
  // @todo Please insert your code here
}


void ElasticWave::ConservedQuantitiesWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i< 3; i++){
    outputQuantities[i] = Q[i];
  }
  for (int i=0; i< 3; i++){
    outputQuantities[i+3] = Q[i+9];
  }
}


