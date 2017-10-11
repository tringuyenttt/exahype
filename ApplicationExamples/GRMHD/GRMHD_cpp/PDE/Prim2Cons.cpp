
#include "PDE.h"
#include <cmath>

using namespace GRMHD;

constexpr bool debug_p2c = false;
#define S(x) printf(#x " = %e\n", x);
#define SI(x) {DFOR(i) S(x(i)); }

void GRMHD::Prim2Cons::perform() {
	// The hydro + magneto exact known prim2cons
	Dens = rho * W;
	DFOR(i) Si.lo(i) = Dens*enth*W*vel.lo(i) + BmagBmag*vel.lo(i) - BmagVel*Bmag.lo(i);
	tau = Dens*enth*W - press + 0.5*(BmagBmag*(1+VelVel) - BmagVel*BmagVel);
}

void GRMHD::Prim2Cons::prepare() {
	constexpr double gamma = GRMHD::Parameters::gamma;
	double epsilon;

	// we compute v_i and B_i on our own. Of course they could in principle also
	// be given by some ID code or so.
	vel.lo   = 0; DFOR(i) CONTRACT(j)  vel.lo(i) +=  vel.up(j) * gam.lo(i,j);
	Bmag.lo  = 0; DFOR(i) CONTRACT(j) Bmag.lo(i) += Bmag.up(j) * gam.lo(i,j);
	
	VelVel   = 0; CONTRACT(i)   VelVel +=  vel.lo(i) *  vel.up(i); // v^2 = v_i * v^i
	BmagBmag = 0; CONTRACT(k) BmagBmag += Bmag.lo(k) * Bmag.up(k); // B^2 = B^j * B_j
	BmagVel  = 0; CONTRACT(k)  BmagVel += Bmag.lo(k) *  vel.up(k); // B^j * v_j
	
	W = 1. / std::sqrt(1.0 - VelVel); // Lorentz factor
	if(debug_p2c) { printf("P2C: "); S(W); S(VelVel); }

	// If we store the pressure (as we do currently in ExaHyPE) in the primitive variables, i.e.
	// V=(rho,vel,press) instead of V=(rho,vel,eps) with internal energy eps, we need to recover
	// eps = eps(press) as the equation of state:
	epsilon = press / (rho * (gamma-1));
	// The specific enthalpy including restmass is now independent of the EOS:
	enth = 1 + epsilon + press / rho;
}

void GRMHD::Prim2Cons::copyFullStateVector() {
	// 1) Assume that Conserved Hydro vars have been set
	// 2) Copy Magneto variables
	copy_magneto(Q);
	copy_adm(Q);
}
