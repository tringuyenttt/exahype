#include "GRMHDSolver_ADERDG.h"
#include <algorithm> // fill_n
#include <cstring> // memset

#include "kernels/KernelUtils.h" // matrix indexing
#include "kernels/GaussLegendreQuadrature.h"

#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
//namespace GRMHD { constexpr int nVar = GRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
#include "PDE/PDE.h"
#include "InitialData/InitialData.h"

#include "GRMHDSolver_ADERDG_Variables.h"
#include "DebuggingHelpers.h"


constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
constexpr int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
constexpr int basisSize = order + 1;
constexpr int nDim = DIMENSIONS;

tarch::logging::Log GRMHD::GRMHDSolver_ADERDG::_log( "GRMHD::GRMHDSolver_ADERDG" );

void GRMHD::GRMHDSolver_ADERDG::init(std::vector<std::string>& cmdlineargs) { // ,  exahype::Parser::ParserView constants) {
	// feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	
	// Initialize initial data
	InitialDataCode::getInstance();
	
}

void GRMHD::GRMHDSolver_ADERDG::adjustPointSolution(const double* const x,const double t,const double dt,double* Q) {
  if (tarch::la::equals(t,0.0)) {
    InitialData(x,t,Q);
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
	
	// as TDIM is 3 and DIMENSIONS is 2, come up with this dirty wrapper:
	double *FT[3], Fz[nVar];
	FT[0] = F[0];
	FT[1] = F[1];
	FT[2] = Fz;
	
	GRMHD::Fluxes(FT, Q).zeroMaterialFluxes();
	//DFOR(d) { zero2Din3D(F[d]); zeroHelpers(F[d]); }
	zeroHelpers(F[0]); zeroHelpers(F[1]);
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

		InitialData(x, ti, Qgp);
		flux(Qgp, F);

		for(int m=0; m < nVar; m++) {
			stateOut[m] += weight * Qgp[m];
			fluxOut[m] += weight * Fs[d][m];
		}
	}
	
	//NVARS(i) printf("stateOut[%d]=%e\n", i, stateOut[i]);
	//NVARS(i) printf("fluxOut[%d]=%e\n", i, fluxOut[i]);
	//std::abort();
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
	//zero2Din3D(S_);
}

void GRMHD::GRMHDSolver_ADERDG::nonConservativeProduct(const double* const Q,const double* const gradQ, double* BgradQ) {
	PDE::NCP ncp(BgradQ);
	
	// as we only have DIMENSIONS == 2 but TDIM == 3, construct a zero gradient
	double gradZ[nVar] = {0.};
	
	// deconstruct the gradient matrix because we use more than the 19 variables of GRMHD.
	const Gradients g(gradQ+0, gradQ+nVar, gradZ);
	PDE pde(Q);
	pde.nonConservativeProduct(g, ncp);
	ncp.zero_adm();
	zeroHelpers(BgradQ);
	//zero2Din3D(BgradQ);
	
	// the NCP must be zero in flat space
	//constexpr double eps = 1e-10;
	//NVARS(i) { if(BgradQ[i]>eps) { printf("BgradQ[%d] = %e\n", i, BgradQ[i]); std::abort(); } }
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
