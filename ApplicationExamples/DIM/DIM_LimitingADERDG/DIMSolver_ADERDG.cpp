#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

#include "PDE.h"

#include "InitialData.h"
// Used for the rieman-solver part
#include "peano/utils/Loop.h"
#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"
tarch::logging::Log DIM::DIMSolver_ADERDG::_log( "DIM::DIMSolver_ADERDG" );


void DIM::DIMSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

void DIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 14 + 0
  // @todo Please implement/augment if required
  if (tarch::la::equals(t,0.0)) {
    // Fortran
    initialdata_(x, &t, Q);
  }
}

void DIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 14 + 0

  // @todo Please implement/augment if required

  const int nVar = DIM::AbstractDIMSolver_ADERDG::NumberOfVariables;
  const int order = DIM::AbstractDIMSolver_ADERDG::Order;
  const int basisSize = order + 1;
  const int nDim = DIMENSIONS;

  double Qgp[nVar], F[nDim][nVar];

  std::memset(stateOut, 0, nVar * sizeof(double));
  std::memset(fluxOut, 0, nVar * sizeof(double));
  
  for(int i=0; i < basisSize; i++)  { // i == time
     const double weight = kernels::gaussLegendreWeights[order][i];
     const double xi = kernels::gaussLegendreNodes[order][i];
     double ti = t + xi * dt;

     initialdata_(x, &ti, Qgp);
     pdeflux_(F[0], F[1], (nDim==3)?F[2]:nullptr, Qgp);
     for(int m=0; m < nVar; m++) {
        stateOut[m] += weight * Qgp[m];
        fluxOut[m] += weight * F[normalNonZero][m];
     }
  }
}

exahype::solvers::Solver::RefinementControl DIM::DIMSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

//*****************************************************************************
//******************************** PDE ****************************************
// To use other PDE terms, specify them in the specification file, delete this 
// file and its header and rerun the toolkit
//*****************************************************************************


void DIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void DIM::DIMSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(
    double* observables,const int numberOfObservables,
    const double* const Q) const {
  assertion(numberOfObservables==1);
  observables[0]=Q[12]; //extract alpha
}

bool DIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
  const double* const solution,
  const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
  const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
  const double t, const double dt) const {
  
//  if ((center[0]-0.5)*(center[0]-0.5)+(center[1]-0.5)*(center[1]-0.5)<0.25*dx[0]*dx[0]) return false;

  if (observablesMin[0] <= 0.0) return false;
  return true;
}


void  DIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


