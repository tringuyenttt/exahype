#include "DummySolver_ADERDG.h"

#include "DummySolver_ADERDG_Variables.h"

#include "Shared.h"

tarch::logging::Log Dummy::DummySolver_ADERDG::_log( "Dummy::DummySolver_ADERDG" );


void Dummy::DummySolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void Dummy::DummySolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t,0.0)) initialdata(x,Q);
}

void Dummy::DummySolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  initialdata(x,stateOut);
  set(fluxOut, 0);
}

exahype::solvers::Solver::RefinementControl Dummy::DummySolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void Dummy::DummySolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  set(lambda, 1.0);
}


void Dummy::DummySolver_ADERDG::flux(const double* const Q,double** F) {
  for(int i=0; i<DIMENSIONS; i++) set(F[i], 0);
}





