#include "MyElasticSolver.h"

#include <memory>

Elastic::MyElasticSolver::MyElasticSolver(int nodesPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):
  exahype::solvers::ADERDGSolver("MyElasticSolver", 12, 0, nodesPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {
  // @todo Please implement/augment if required
}



void Elastic::MyElasticSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)

  double* f = F[0];
  double* g = F[1];
  double* h = F[2];

   double rho;       /* density */
  double lam;       /* 1st Lame parameter */
  double mu;        /* 2nd Lame parameter */
  
  rho = Q[9];
  lam = Q[10];
  mu = Q[11];
  
  // @todo Please implement
  // f
  f[0] = 1.0/rho*Q[3];
  f[1] = 1.0/rho*Q[6];
  f[2] = 1.0/rho*Q[7];
  f[3] = (2.0*mu + lam)*Q[0];
  f[4] = lam*Q[0];
  f[5] = lam*Q[0];
  f[6] = mu*Q[1];
  f[7] = mu*Q[2];
  f[8] = 0.0;
  f[9] = 0.0;
  f[10] = 0.0;
  f[11] = 0.0;
  // g
  // @todo Please implement
  g[0] = 1.0/rho*Q[6];
  g[1] = 1.0/rho*Q[4];
  g[2] = 1.0/rho*Q[8];
  g[3] = lam*Q[1];
  g[4] = (2.0*mu + lam)*Q[1];
  g[5] = lam*Q[1];
  g[6] = mu*Q[0];
  g[7] = 0.0;
  g[8] = mu*Q[2];
  g[9] = 0.0;
  g[10] = 0.0;
  g[11] = 0.0;
  // h
  // @todo Please implement
  h[0] = 1.0/rho*Q[7];
  h[1] = 1.0/rho*Q[8];
  h[2] = 1.0/rho*Q[5];
  h[3] = lam*Q[2];
  h[4] = lam*Q[2];
  h[5] = (2.0*mu + lam)*Q[2];
  h[6] = 0.0;
  h[7] = mu*Q[0];
  h[8] = mu*Q[1];
  h[ 9] = 0.0;
  h[10] = 0.0;
  h[11] = 0.0;
}



void Elastic::MyElasticSolver::source(const double* const Q, double* S) {
  // Number of variables = 12 + 0
  // @todo Please implement
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = 0.0;
  S[6] = 0.0;
  S[7] = 0.0;
  S[8] = 0.0;
  S[9] = 0.0;
  S[10] = 0.0;
  S[11] = 0.0;
}



void Elastic::MyElasticSolver::boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)

  // @todo Please implement
  // fluxOut
  fluxOut[ 0] = fluxIn[ 0];
  fluxOut[ 1] = fluxIn[ 1];
  fluxOut[ 2] = fluxIn[ 2];
  fluxOut[ 3] = 0*fluxIn[ 3];
  fluxOut[ 4] = 0*fluxIn[ 4];
  fluxOut[ 5] = 0*fluxIn[ 5];
  fluxOut[ 6] = 0*fluxIn[ 6];
  fluxOut[ 7] = 0*fluxIn[ 7];
  fluxOut[ 8] = fluxIn[ 8];
  fluxOut[ 9] = fluxIn[ 9];
  fluxOut[10] = fluxIn[10];
  fluxOut[11] = fluxIn[11];
  // stateOut
  // @todo Please implement
  stateOut[ 0] = 0*stateIn[ 0];
  stateOut[ 1] = 0*stateIn[ 1];
  stateOut[ 2] = 0*stateIn[ 2];
  stateOut[ 3] = stateIn[ 3];
  stateOut[ 4] = stateIn[ 4];
  stateOut[ 5] = stateIn[ 5];
  stateOut[ 6] = stateIn[ 6];
  stateOut[ 7] = stateIn[ 7];
  stateOut[ 8] = stateIn[ 8];
  stateOut[ 9] = stateIn[ 9];
  stateOut[10] = stateIn[10];
  stateOut[11] = stateIn[11];

  if (faceIndex == 0){
    fluxOut[ 0] = fluxIn[ 0];
    fluxOut[ 1] = fluxIn[ 1];
    fluxOut[ 2] = fluxIn[ 2];
    fluxOut[ 3] = fluxIn[ 3];
    fluxOut[ 4] = fluxIn[ 4];
    fluxOut[ 5] = fluxIn[ 5];
    fluxOut[ 6] = fluxIn[ 6];
    fluxOut[ 7] = fluxIn[ 7];
    fluxOut[ 8] = fluxIn[ 8];
    fluxOut[ 9] = fluxIn[ 9];
    fluxOut[10] = fluxIn[10];
    fluxOut[11] = fluxIn[11];
    // stateOut
    // @todo Please implement
    stateOut[ 0] = stateIn[ 0];
    stateOut[ 1] = stateIn[ 1];
    stateOut[ 2] = stateIn[ 2];
    stateOut[ 3] = stateIn[ 3];
    stateOut[ 4] = stateIn[ 4];
    stateOut[ 5] = stateIn[ 5];
    stateOut[ 6] = stateIn[ 6];
    stateOut[ 7] = stateIn[ 7];
    stateOut[ 8] = stateIn[ 8];
    stateOut[ 9] = stateIn[ 9];
    stateOut[10] = stateIn[10];
    stateOut[11] = stateIn[11];
    
    fluxOut[0] = 0.0*fluxIn[0];
    fluxOut[1] = 0.0*fluxIn[1];
    fluxOut[2] = 0.0*fluxIn[2];
    
 
    stateOut[3] = 0.0*stateIn[3];
    stateOut[6] = 0.0*stateIn[6];
    stateOut[7] = 0.0*stateIn[7];
    
  }
}



void Elastic::MyElasticSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  // @todo Please implement

   double rho;        /* density */
  double lam;        /*  1st Lame parameter */
  double mu;         /* 2nd Lame parameter */
  
  rho = Q[9];
  lam = Q[10];
  mu = Q[11];

  // extract p-wave speed (cp) and s-wave speed (cs)
  double cp = std::sqrt((2.0*mu + lam)/rho);
  double cs = std::sqrt(mu/rho);
  
  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = -cs;
  lambda[3] = cp;
  lambda[4] = cs;
  lambda[5] = cs;
  lambda[6] = 0.0;
  lambda[7] = 0.0;
  lambda[ 8] = 0.0;
  lambda[ 9] = 0.0;
  lambda[10] = 0.0;
  lambda[11] = 0.0;
  
}



bool Elastic::MyElasticSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
  // @todo Please implement
   if (tarch::la::equals(t, 0.0)) {
    return true;
  }
  return false;
}



void Elastic::MyElasticSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 12 (#unknowns + #parameters)
  // @todo Please implement
   double rho;        /* density */
  double cp;        /*  p-wave speed */
  double cs;         /* s-wave speed */

  rho = 2.6;      // gm/cm^3
  cs = 2.0;       // km/s
  cp = 4.0;       // km/s
  
  double x0 = 5;
  double y0 = 10;
  double z0 = 10;
  
  Q[0] = std::exp(-((x[0]-x0)*(x[0]-x0) + (x[1]-y0)*(x[1]-y0) + (x[2]-z0)*(x[2]-z0) ) /
		  (1.0 * 1)) *10.0;
  Q[1] = std::exp(-((x[0]-x0)*(x[0]-x0) + (x[1]-y0)*(x[1]-y0) + (x[2]-z0)*(x[2]-z0) ) /
		  (1 * 1)) *10.0;
  Q[2] = std::exp(-((x[0]-x0)*(x[0]-x0) + (x[1]-y0)*(x[1]-y0) + (x[2]-z0)*(x[2]-z0) ) /
		  (1 * 1)) *10.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
  Q[5] = 0.0;
  Q[6] = 0.0;
  Q[7] = 0.0;
  Q[8] = 0.0;
  
  Q[9] = rho; 
  Q[10] = rho*(cp*cp - 2*cs*cs) ;  // lambda (MPa)
  Q[11] = rho*cs*cs;               // mu (MPa)

  if (x[0] > 1.0){  
    rho = 2.7;      // gm/cm^3
    cs = 3.464;       // km/s
    cp = 6.0;       // km/s

    Q[9] = rho;                     // rho 
    Q[10] = rho*(cp*cp - 2*cs*cs) ;  // lambda (MPa)
    Q[11] = rho*cs*cs;               // mu     (MPa) 
  }
  
}



exahype::solvers::Solver::RefinementControl Elastic::MyElasticSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
  // @todo Please implement
  return exahype::solvers::Solver::RefinementControl::Keep;
}


