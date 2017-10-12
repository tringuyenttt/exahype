#ifndef __GRMHD_FOR_EXAHYPE_ADAPTER__
#define __GRMHD_FOR_EXAHYPE_ADAPTER__


/**
 * Coupling of the SVEC GRMHD PDE to ExaHyPE.
 **/

#include "PDE/PDE.h"

namespace SVEC {
namespace GRMHDAdapterForExaHyPE {

#ifdef TEST_NEW_PDE_AUTONOMOUSLY
	// in order to autonomously test/copmile this C++ file:
	//#define TDIM 3
	constexpr int nVar = 19;
	//int main() { return 0; }
#else
	// If you include to ExaHyPE instead:
	#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
	#include "AbstractGRMHDSolver_FV.h" // Defines: nVar
	constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables;  // ensure this is 19 or so
	//namespace SVEC { constexpr int nVar = 19; }
#endif

//namespace SVEC { constexpr int nVar = 42; }
#define NVARS(i) for(int i=0; i < SVEC::GRMHDAdapterForExaHyPE::nVar; i++)

	/********* CONSERVATIVE PART ************/
	
	// Represents fluxes after being computes
	struct Fluxes {
		GRMHD::Flux F[DIMENSIONS];
		#if DIMENSIONS == 2
		Fluxes(double *fx, double *fy) : F{fx,fy} {}
		Fluxes(double** f) : Fluxes(f[0],f[1]) {}
		Fluxes(double *fx, double *fy, double *fz) : Fluxes(fx,fy) {} ///< for convention
		#elif DIMENSIONS == 3
		Fluxes(double *fx, double *fy) : F{fx,fy} {}
		Fluxes(double** f) : Fluxes(f[0],f[1]) {}
		#endif
		GRMHD::Flux& operator[](int i) { return F[i]; }
		void zeroMaterialFluxes() {
			for(int d=0; d<DIMENSIONS; d++) F[d].zero_adm();
		}
	};

	Fluxes flux(const double* const Q, double** F) {
		GRMHD::PDE::Conserved pde(Q);
		Fluxes f(F);
		for(int d=0; d<DIMENSIONS; d++) pde.flux(f[d], d);
		return f;
	}
	
	// we also can work with splitted fluxes:
	Fluxes flux(const double* const Q, double *Fx, double *Fy, double *Fz) {
		GRMHD::PDE::Conserved pde(Q);
		Fluxes f(Fx,Fy,Fz);
		for(int d=0; d<DIMENSIONS; d++) pde.flux(f[d], d);
		return f;
	}
	
	// Gradients coming from double pointers
	struct Gradients {
		const GRMHD::Gradient dir[TDIM];
		// Access an element in some direction
		const GRMHD::Gradient& operator[](int i) const { return dir[i]; }
		// Access an element, indicating that it's partial_i, not partial^i
		const GRMHD::Gradient& lo(int i) const { return dir[i]; }
		
		// TODO IMPORTANT:
		// Consider fact that DIMENSIONS == 2 but TDIM ==3.
		
		#if TDIM == 3
		/// individual vectors
		Gradients(const double* const Qx, const double* const Qy, const double* const Qz) : dir{Qx,Qy,Qz} {}
		/// similar form as the fluxes (matrix)
		Gradients(const double** const gradQ) : Gradients(gradQ[0],gradQ[1],gradQ[2]) {}
		/// linear array
		Gradients(const double* const gradQ) : Gradients(gradQ+0,gradQ+nVar,gradQ+2*nVar) {}
		#elif TDIM == 2
		/// individual vectors
		Gradients(const double* const Qx, const double* const Qy) : dir{Qx,Qy} {}
		/// Similar form as the fluxes
		Gradients(const double** const gradQ) : Gradients(gradQ[0],gradQ[1]) {}
		/// Linear array
		Gradients(const double* const gradQ) : Gradients(gradQ+0,gradQ+nVar) {}
		/// For convenience a similar interface as what TDIM==3 provides:
		Gradients(const double* const Qx, const double* const Qy, const double* const Qz) : Gradients(Qx,Qy) {}
		#endif
	};
	
	/********* NONCONSERVATIVE PART ************/

	GRMHD::NCP nonConservativeProduct(const double* const Q, const double* const gradQ, double* BgradQ) {
		GRMHD::PDE pde(Q);
		GRMHD::NCP n(BgradQ);
		Gradients g(gradQ);
		pde.nonConservativeProduct<Gradients>(g, n);
		return n;
	}
	
	GRMHD::NCP nonConservativeProduct(const double* const Q, const double* const Qx, const double* const Qy, const double* const Qz, double* BgradQ) {
		GRMHD::PDE pde(Q);
		GRMHD::NCP n(BgradQ);
		Gradients g(Qx,Qy,Qz);
		pde.nonConservativeProduct<Gradients>(g, n);
		return n;
	}
	
	GRMHD::Source algebraicSource(const double* const Q, double* S) {
		GRMHD::PDE pde(Q);
		GRMHD::Source s(S);
		pde.algebraicSource(s);
		return s;
	}
	
	GRMHD::Source fusedSource(const double* const Q, const double* const gradQ, double* S) {
		GRMHD::PDE pde(Q);
		GRMHD::Source s(S);
		Gradients g(gradQ);
		pde.fusedSource<Gradients>(g, s);
		return s;
	}

	GRMHD::Source fusedSource(const double* const Q, const double* const Qx, const double* const Qy, const double* const Qz, double* S) {
		GRMHD::PDE pde(Q);
		GRMHD::Source s(S);
		Gradients g(Qx,Qy,Qz);
		pde.fusedSource<Gradients>(g, s);
		return s;
	}
	
} // ns GRMHDAdapterForExaHyPE
} // ns SVEC

#endif /* __GRMHD_FOR_EXAHYPE_ADAPTER__ */
