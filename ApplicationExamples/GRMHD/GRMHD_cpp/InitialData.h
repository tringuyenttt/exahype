
#include "PDE/PDE.h"
#include <cmath> 

using namespace GRMHD;

struct InitialState : public GRMHDSystem::Shadow, public Hydro::Primitives::Stored {
	InitialState(double* Q) : GRMHDSystem::Shadow(Q), Hydro::Primitives::Stored() {
		// by default, set NaNs for uninitialized values to catch mistakes when setting ID
		Dens = tau = phi = rho = press = alpha = NAN;
		DFOR(i) { Bmag.up(i) = Si.lo(i) = vel.up(i) = beta.up(i) = NAN; }
		DFOR(i) DFOR(j) gam.lo(i,j) = NAN;
	}
};

/// Defaults to vacuum
struct VacuumInitialData : public InitialState {
	static constexpr double rho0 = 1.;
	static constexpr double p0   = 1.;
	
	void flatSpace() {
		// ADM base: Flat space
		alpha = 1.0;
		DFOR(i) beta.up(i) = 0;
		SYMFOR(i,j) gam.lo(i,j) = sym::delta(i,j);
	}
	
	VacuumInitialData(double* Q) : InitialState(Q) {
		// zero first all variables. This is useful for a 2D workaround as 3D variables
		// are stored but not initialized (only 2D storage is actually used)
		// NVARS(i) Q[i] = 0;
		
		// For test, do not set anything
		//return;
		
		// c2p invariants:
		DFOR(i) Bmag.up(i) = 0;
		phi = 0; // damping term
		
		// primitives:
		rho = rho0;
		press = p0;
		DFOR(i) vel.up(i) = 0;
		
		// ADM
		flatSpace();
		
		// Don't forget to call the Cons2Prim afterwards.
	}
};

struct AlfenWave : public VacuumInitialData {
	AlfenWave(const double* const x, const double t, double* Q);
};

struct AlfenWaveCons : public AlfenWave {
	double V[100];
	AlfenWaveCons(const double* const x, const double t, double* Q) : AlfenWave(x,t,V) {
			GRMHD::Prim2Cons(Q, V).copyFullStateVector();
		}
};
