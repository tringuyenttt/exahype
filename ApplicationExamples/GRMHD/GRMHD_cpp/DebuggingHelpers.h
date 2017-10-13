// shared between ADERDG and FV solvers

#include "GRMHDSolver_ADERDG_Variables.h"
#include <fenv.h> // enable nan tracker

constexpr double magicCheck = 123456789.0123456789;

inline void enableNanCatcher() {
	feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
}

inline void zeroHelpers(double* Q) {
	// debugging variables
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts Qi;
	Q[Qi.sqdetg] = 0;
	for(int i=0;i<2;i++) Q[Qi.spare+i] = 0;
	var.check() = 0;
}

inline void zero2Din3D(double* Q) {
	// don't do that any more.
	return;
	/*
	
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
	*/
}

// The shared ID, moved here temporarily
/*
inline void InitialData(const double* const x,const double t,double* Q) {
  // Number of variables    = 23 + #parameters
  NVARS(i) Q[i] = NAN; // to find problems

  // currently, the C++ AlfenWave spills out primitive data
  //  double V[nVar];
  //  AlfenWave id(x,t,V);
  //  GRMHD::Prim2Cons(Q, V).copyFullStateVector();
  AlfenWaveCons(x,t,Q);

  // also store the positions for debugging
  GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
  GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts Qi;
  for(int i=0; i<DIMENSIONS; i++) var.pos(i) = x[i];
  if(DIMENSIONS<3) Q[Qi.pos+2] = -1;
  var.check() = magicCheck;
  //zero2Din3D(Q);

  NVARS(i) { if(!std::isfinite(Q[i])) { printf("Qid[%d] = %e\n", i, Q[i]); std::abort(); } }
}
*/
  
// Detection of unphysical states. In these cases, the user PDE functions shall never be called.
// We workaround by returning some kind of "neutral" values which go well with the scheme.

inline bool isAllZero(const double* const Q) {
	bool allzero=true; NVARS(i) { if(Q[i]!=0) allzero=false; }
	return allzero;
}

// Check whether Q holds the vacuum spacetime up to some uncertainty
inline bool holdsVacuumSpacetime(const double* const Q) {
	using namespace SVEC;
	using namespace tarch::la;
	
	ADMBase::ConstShadow adm(Q);
	bool correct = equals(adm.alpha, 1);
	DFOR(i) correct &= equals(adm.beta.up(i), 0);
	CONTRACT2(i,j) correct &= equals(adm.gam.lo(i,j), delta(i,j));
	return correct;
}


inline void fail(const std::string msg) {
	fputs((msg + "\n").c_str(), stderr);
	std::abort();
}

inline bool isUnphysical(const double* const Q) {
	bool isUnphysical = isAllZero(Q);// || !holdsVacuumSpacetime(Q);
	if(isUnphysical) fail("is unphysical");
}



/*
inline bool metrikIsUseful() {
	PDE pde(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::ReadOnlyVariables var(Q);
	constexpr double eps = 1e-8;
	// ExaHyPE workaround:
	// if the input has zeros at weird places, don't do anything
	return !(std::abs(pde.gam.det) < eps || var.check()!=magicCheck);
}
*/
