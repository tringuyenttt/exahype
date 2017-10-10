#include "GRMHDSolver_ADERDG.h"
#include <algorithm> // fill_n
#include <cstring> // memset
#include <fenv.h> // enable nan tracker

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
//namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
#include "PDE/PDE.h"
#include "InitialData.h"

#include "GRMHDSolver_ADERDG_Variables.h"


constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log( "GRMHD::GRMHDSolver_ADERDG" );

constexpr double magicCheck = 123456789;

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) { // ,  exahype::Parser::ParserView constants) {

	// feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
}

void zeroHelpers(double* Q) {
	// debugging variables
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	DFOR(i) var.pos(i) = 0;
	var.check() = 0;
}

void zero2Din3D(double* Q) {
	// 3D components which should not be in 2D
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts pos;
	
	// 3rd component of vectors
	Q[pos.vel+2] = 0;
	Q[pos.B  +2] = 0;
	Q[pos.shift+2] = 0;
	// 3d-components of tensor
	for(int i=0; i<3; i++) Q[pos.gij+tensish::sym::index(2,i)] = 0;
	 // 3rd component of debugging coordinate vector
	Q[pos.pos+2] = 0;
}

void /*GRMHD::GRMHDSolver_ADERDG:: */ initialData(const double* const x,const double t,const double dt,double* Q) {
  // Number of variables    = 23 + #parameters
  NVARS(i) Q[i] = NAN; // to find problems

  // currently, the C++ AlfenWave spills out primitive data
  /*
    double V[nVar];
    AlfenWave id(x,t,V);
    GRMHD::Prim2Cons(Q, V).copyFullStateVector();
   */
  AlfenWaveCons(x,t,Q);

  // also store the positions for debugging
  GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
  DFOR(i) var.pos(i) = x[i];
  var.check() = magicCheck;
  zero2Din3D(Q);

  NVARS(i) { if(!std::isfinite(Q[i])) { printf("Qid[%d] = %e\n", i, Q[i]); std::abort(); } }
}

void GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t,0.0)) {
    initialData(x,t,dt,Q);
  }
}

void GRMHD::GRMHDSolver_ADERDG::eigenvalues(const double* const Q,const int d,double* lambda) {
	// Provide NVARS eigenvalues
	//PDE::eigenvalues(Q, d, lambda);
	
	// Just set all eigenvalues to zero.
	vec::shadow<nVar> Lambda(lambda);
	Lambda = 1.0;
}


void GRMHD::GRMHDSolver_ADERDG::flux(const double* const Q,double** F) {
	GRMHD::Fluxes(F, Q).zeroMaterialFluxes();
	DFOR(d) { zero2Din3D(F[d]); zeroHelpers(F[d]); }
}


void GRMHD::GRMHDSolver_ADERDG::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int d,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
	// for debugging, to make sure BC are set correctly
	/*
	double snan = std::numeric_limits<double>::signaling_NaN();
	double weird_number = -1.234567;
	std::memset(stateOut, weird_number, nVar * sizeof(double));
	std::memset(fluxOut,  weird_number, nVar * sizeof(double));
	*/
	
	// employ time-integrated exact BC for AlfenWave.

	double Qgp[nVar], Fs[nDim][nVar], *F[nDim];
	for(int dd=0; dd<nDim; dd++) F[dd] = Fs[dd];
	// zeroise stateOut, fluxOut
	for(int m=0; m<nVar; m++) {
		stateOut[m] = 0;
		fluxOut[m] = 0;
	}
	for(int i=0; i < basisSize; i++)  { // i == time
		const double weight = kernels::gaussLegendreWeights[order][i];
		const double xi = kernels::gaussLegendreNodes[order][i];
		double ti = t + xi * dt;

		initialData(x, ti, dt, Qgp);
		flux(Qgp, F);

		for(int m=0; m < nVar; m++) {
			stateOut[m] += weight * Qgp[m];
			fluxOut[m] += weight * Fs[d][m];
		}
	}
}


exahype::solvers::Solver::RefinementControl GRMHD::GRMHDSolver_ADERDG::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}

void GRMHD::GRMHDSolver_ADERDG::algebraicSource(const double* const Q,double* S_) {
	PDE::Source S(S_);
	PDE(Q).algebraicSource(S);
	S.zero_adm();
	zeroHelpers(S_);
	zero2Din3D(S_);
}

void GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ, double* BgradQ) {
	PDE::NCP ncp(BgradQ);
	
	// deconstruct the gradient matrix because we use more than the 19 variables of GRMHD.
	const Gradients g(gradQ+0, gradQ+nVar, gradQ+2*nVar);
	PDE pde(Q);
	pde.nonConservativeProduct(g, ncp);
	ncp.zero_adm();
	zeroHelpers(BgradQ);
	zero2Din3D(BgradQ);
	
	// the NCP must be zero in flat space
	constexpr double eps = 1e-10;
	NVARS(i) { if(BgradQ[i]>eps) { printf("BgradQ[%d] = %e\n", i, BgradQ[i]); std::abort(); } }
}

// The optimized fusedSource
/*
void GRMHD::GRMHDSolver_ADERDG::fusedSource(const double* const Q, const double* const gradQ, double* S_) {
	PDE::Source S(S_);
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::ReadOnlyVariables var(Q);
	constexpr double eps = 1e-8;
	// ExaHyPE workaround:
	// if the input has zeros at weird places, don't do anything
	if(std::abs(pde.gam.det) < eps || var.check()!=magicCheck) {
		printf("Weird FusedSource input (det=%e), skipping. ",pde.gam.det);
		//SYMFOR(i,j) printf("gam(%d,%d)=%e\n", i,j, pde.gam.lo(i,j));
		std::cout << "pos: " << var.pos() << " check: " << var.check() <<  std::endl;
		// set everything to NANs
		NVARS(m) S_[m] = NAN;
		return;
	} else {
		//printf("Good FusedSource input (det=%e)\n",pde.gam.det);
		//NVARS(m) printf("FusedSource input: Q[%d]=%e\n", m, Q[m]);
		
		
		//std::cout << "Good FusedSource, x= " << var.pos() <<  std::endl;
	}
	pde.RightHandSide(gradQ, S);
	S.zero_adm();
	zeroHelpers(S_);
	zero2Din3D(S_);
}
*/
