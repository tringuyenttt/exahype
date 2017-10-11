#include "DIMSolver_ADERDG.h"

#include "DIMSolver_ADERDG_Variables.h"

#include "PDE.h"

#include "InitialData.h"

#include "peano/utils/Loop.h"

// #include "Limiter.h"

#include <cstring> // memset
#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

tarch::logging::Log DIM::DIMSolver_ADERDG::_log( "DIM::DIMSolver_ADERDG" );

bool DIM::DIMSolver_ADERDG::isPhysicallyAdmissible(
  const double* const solution,
  const double* const observablesMin,const double* const observablesMax,const int numberOfObservables,
  const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx,
  const double t, const double dt) const {
  // True == NoLimiter, False == Limiter
  const ReadOnlyVariables qmin(observablesMin);

// Limiter based on the center of the cell
const double rad = 0.25;

// points
int outsidePoints = 0;

tarch::la::Vector<DIMENSIONS,double> offset = center-0.5*dx;
dfor2(p)
  tarch::la::Vector<DIMENSIONS,double> corner=offset;
  #if DIMENSIONS==3
  corner[2]+=p[2]*dx[2];
  #endif
  corner[1]+=p[1]*dx[1];
  corner[0]+=p[0]*dx[0];
  
  outsidePoints += (tarch::la::norm2(corner)>rad) ? 1 : 0;
enddforx

if (outsidePoints>0 && outsidePoints<TWO_POWER_D) 
{
  return false;
}
else
{
  return true; 
}

//if(tarch::la::norm2(center)<0.32 &&
//tarch::la::norm2(center)>0.2 )
//{
//  return false;
//}
//else
//{
//  return true;
//}

/*
// Physical admissibility criteria done in C++
if(QMin[ 12]<0.99999 && QMax[ 12]>0.00001)
{
  return false;
}else{
  return true; 
}
*/
/*
// Physical admissibility critedia done in fortran
int status;

PDEAssurePositivity_(&status, QMin, QMax);

if(status<0)
{
  return false;
}else{
  return true; 
}
*/
}


void DIM::DIMSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue DIM::DIMSolver_ADERDG::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
 /*
  Q[ 0] = 0.0;
  Q[ 1] = 0.0;
  Q[ 2] = 0.0;
  Q[ 3] = 0.0;
  Q[ 4] = 0.0;
  Q[ 5] = 0.0;
  Q[ 6] = 0.0;
  Q[ 7] = 0.0;
  Q[ 8] = 0.0;
  Q[ 9] = 0.0;
  Q[10] = 0.0;
  Q[11] = 0.0;
  Q[12] = 0.0;
  Q[13] = 0.0;
*/
  // Fortran
  if (tarch::la::equals(t,0.0)) {
    initialdata_(x, &t, Q);
  }
}

void  __attribute__((optimize("O0"))) DIM::DIMSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters
  
  // @todo Please implement/augment if required
/*
  lambda[ 0] = 1.0;
  lambda[ 1] = 1.0;
  lambda[ 2] = 1.0;
  lambda[ 3] = 1.0;
  lambda[ 4] = 1.0;
  lambda[ 5] = 1.0;
  lambda[ 6] = 1.0;
  lambda[ 7] = 1.0;
  lambda[ 8] = 1.0;
  lambda[ 9] = 1.0;
  lambda[10] = 1.0;
  lambda[11] = 1.0;
  lambda[12] = 1.0;
  lambda[13] = 1.0;
*/
  double nv[3] = {0.};
  nv[d] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}


void __attribute__((optimize("O0"))) DIM::DIMSolver_ADERDG::flux(const double* const Q,double** F) {
  F[0][ 0] = 0.0;
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;
  F[0][ 5] = 0.0;
  F[0][ 6] = 0.0;
  F[0][ 7] = 0.0;
  F[0][ 8] = 0.0;
  F[0][ 9] = 0.0;
  F[0][10] = 0.0;
  F[0][11] = 0.0;
  F[0][12] = 0.0;
  F[0][13] = 0.0;

  F[1][ 0] = 0.0;
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
  F[1][ 5] = 0.0;
  F[1][ 6] = 0.0;
  F[1][ 7] = 0.0;
  F[1][ 8] = 0.0;
  F[1][ 9] = 0.0;
  F[1][10] = 0.0;
  F[1][11] = 0.0;
  F[1][12] = 0.0;
  F[1][13] = 0.0;
}

void DIM::DIMSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 14 + #parameters

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
        fluxOut[m] += weight * F[d][m];
     }
  }
}



exahype::solvers::Solver::RefinementControl DIM::DIMSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
/*
  if( (center[0]*center[0]+center[1]*center[1]) < 0.25*0.25 )
  {
    if(dx[0]>getMaximumMeshSize()/9)
    {
      return exahype::solvers::Solver::RefinementControl::Refine;
    }
  }
  return exahype::solvers::Solver::RefinementControl::Keep;
*/
}


void  __attribute__((optimize("O0"))) DIM::DIMSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {

  pdencp_(BgradQ, Q, gradQ);


/*
  BgradQ[ 0] = 0.0;
  BgradQ[ 1] = 0.0;
  BgradQ[ 2] = 0.0;
  BgradQ[ 3] = 0.0;
  BgradQ[ 4] = 0.0;
  BgradQ[ 5] = 0.0;
  BgradQ[ 6] = 0.0;
  BgradQ[ 7] = 0.0;
  BgradQ[ 8] = 0.0;
  BgradQ[ 9] = 0.0;
  BgradQ[10] = 0.0;
  BgradQ[11] = 0.0;
  BgradQ[12] = 0.0;
  BgradQ[13] = 0.0;
*/  
}


void DIM::DIMSolver_ADERDG::coefficientMatrix(const double* const Q,const int d,double* Bn) {

  double nv[3] = {0.};
  nv[d] = 1;
  pdematrixb_(Bn, Q, nv);

/*
  Bn[ 0] = 0.0;
  Bn[ 1] = 0.0;
  Bn[ 2] = 0.0;
  Bn[ 3] = 0.0;
  Bn[ 4] = 0.0;
  Bn[ 5] = 0.0;
  Bn[ 6] = 0.0;
  Bn[ 7] = 0.0;
  Bn[ 8] = 0.0;
  Bn[ 9] = 0.0;
  Bn[10] = 0.0;
  Bn[11] = 0.0;
  Bn[12] = 0.0;
  Bn[13] = 0.0;
*/
}

