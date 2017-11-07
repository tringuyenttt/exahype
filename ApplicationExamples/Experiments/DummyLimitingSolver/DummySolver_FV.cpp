#include "DummySolver_FV.h"

#include "DummySolver_FV_Variables.h"

#include "Shared.h"

tarch::logging::Log Dummy::DummySolver_FV::_log( "Dummy::DummySolver_FV" );


void Dummy::DummySolver_FV::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void Dummy::DummySolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t,0.0)) initialdata(x,Q);
}

void Dummy::DummySolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  set(lambda, 1.0);
}

void Dummy::DummySolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* stateOutside) {
  initialdata(x,stateOutside);
}


void Dummy::DummySolver_FV::flux(const double* const Q,double** F) {
  for(int i=0; i<DIMENSIONS; i++) set(F[i], 0);
}




