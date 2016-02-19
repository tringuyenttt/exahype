#include "MyEulerSolver.h"



Euler3d::MyEulerSolver::MyEulerSolver( int kernelNumber):
  exahype::solvers::Solver("MyEulerSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,3+1) {
  // @todo Please implement/augment if required
}



int Euler3d::MyEulerSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}
void Euler3d::MyEulerSolver::flux(const double* const Q, double* f, double* g, double* h) {

  // Dimensions             = 3
  // Number of variables    = 5

  const double GAMMA = 1.4;
  
  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1)*( Q[4] - 0.5* (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) * irho );

  // f
  // @todo Please implement
  f[0] = Q[1];
  f[1] = irho*Q[1]*Q[1] + p;
  f[2] = irho*Q[1]*Q[2];
  f[3] = irho*Q[1]*Q[3];
  f[4] = irho*Q[1]*(Q[4]+p);
  // g
  // @todo Please implement
  g[0] = Q[2];
  g[1] = irho*Q[2]*Q[1];
  g[2] = irho*Q[2]*Q[2] + p;
  g[3] = irho*Q[2]*Q[3];
  g[4] = irho*Q[2]*(Q[4]+p);
  // h
  // @todo Please implement
  h[0] = Q[3];
  h[1] = irho*Q[3]*Q[1];
  h[2] = irho*Q[3]*Q[2];
  h[3] = irho*Q[3]*Q[3] + p;
  h[4] = irho*Q[3]*(Q[4]+p);
}



void Euler3d::MyEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 5
  // @todo Please implement
  lambda[0] = 0.0;
  lambda[1] = 0.0;
  lambda[2] = 0.0;
  lambda[3] = 0.0;
  lambda[4] = 0.0;
}



void Euler3d::MyEulerSolver::initialValues(const double* const x, double* Q) {
  // Dimensions             = 3
  // Number of variables    = 5
  // @todo Please implement
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}



