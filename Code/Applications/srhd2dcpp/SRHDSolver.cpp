#include "SRHDSolver.h"
/**
 * This is the pure C++ 2D SRHD Solver.
 * (Can of course easily be used for 3D, too)
 **/

#include <math.h>
#include "SRHDApplication.h"


srhd2d::SRHDSolver::SRHDSolver( int kernelNumber):
  exahype::solvers::Solver("SRHDSolver",exahype::solvers::Solver::ADER_DG,kernelNumber,5,3+1,exahype::solvers::Solver::GlobalTimeStepping) {
  // @todo Please implement/augment if required
}



int srhd2d::SRHDSolver::getMinimumTreeDepth() const {
  // @todo Please implement
  return 3;
}



void srhd2d::SRHDSolver::flux(const double* const Q, double* f, double* g) {
	// Dimensions             = 2
	// Number of variables    = 5 (#unknowns + #parameters)
	
	double V[nVar];
	cons2prim(V, Q);
	SHORTHANDS(V);
	
	double v2 = vx*vx + vy*vy + vz*vz;
	double lf = 1.0/sqrt(1.0 - v2);
	double w = rho + eos_gamma*p/(eos_gamma-1.0);
	double ww = w*lf*lf;
	
	f[0] = vx*rho*lf;
	f[1] = ww*vx*vx + p;
	f[2] = ww*vx*vy;
	f[3] = ww*vx*vz;
	f[4] = ww*vx - f[0];

	g[0] = vy*rho*lf;
	g[1] = ww*vy*vx;
	g[2] = ww*vy*vy + p;
	g[3] = ww*vy*vz;
	g[4] = ww*vy - g[0];
	
	/*
	h[0] = vz*rho*lf;
	h[1] = ww*vz*vx;
	h[2] = ww*vz*vy;
	h[3] = ww*vz*vz + p;
	h[4] = ww*vz - h[0];
	*/
}



void srhd2d::SRHDSolver::eigenvalues(const double* const Q, const int nnzi, double* lambda) {
	// nnzi = normalNonZeroIndex
	// Dimensions             = 2
	// Number of variables    = 5 (#unknowns + #parameters)

	double V[nVar];
	cons2prim(V, Q);
	SHORTHANDS(V);

	double w = rho + eos_gamma*p/(eos_gamma-1.0);
	double cs2 = eos_gamma*p/w;
	double v2 = vx*vx + vy*vy + vz*vz;
	double den = 1.0/(1.0 + v2*cs2);
	double vn = Q[nnzi+1];
	// Question: How does this normalNonZeroIndex work, in contrast to the nv in Fortran's PDEEigenvalues(Lambda,Q,nv)?
	// double vnsum = 0; for(int j=0;j<3;j++) vnsum += Q[i+j]*Q[i+j];
	// double u = (vnsum == 0.) ? sqrt(v2) : vn;
	double u = vn;
	
	lambda[0] = ( u*(1.0-cs2)-sqrt( cs2*(1.0-v2)*( (1.0-v2*cs2) - u*u*(1.0-cs2) )) )*den;
	lambda[1] = 0.0;
	lambda[2] = 0.0;
	lambda[3] = 0.0;
	lambda[4] = ( u*(1.0-cs2)+sqrt( cs2*(1.0-v2)*( (1.0-v2*cs2) - u*u*(1.0-cs2) )) )*den;
}



bool srhd2d::SRHDSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {
	return tarch::la::equals(t, 0.0);
}



void srhd2d::SRHDSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
	// Dimensions             = 2
	// Number of variables    = 5 (#unknowns + #parameters)
	// @todo Please implement
	
	// Riemann test
	if (tarch::la::equals(t, 0.0)) {
		double V[nVar];
		SHORTHANDS(V);

		if(x[0] < 0.5) {
			rho = 1.0;
			p = 1.0;
		} else {
			rho = 0.1;
			p = 0.1;
		}
		vx = vy = vz = 0.;

		prim2cons(V, Q);
	}
}



exahype::solvers::Solver::RefinementControl srhd2d::SRHDSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {
	// @todo Please implement
	return exahype::solvers::Solver::Keep;
}


