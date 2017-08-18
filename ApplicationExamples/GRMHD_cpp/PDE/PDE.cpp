/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * Read PDE.h for further documentation.
 * Written at 2017-08-06 by SvenK.
 **/

#ifdef TEST_NEW_PDE_AUTONOMOUSLY
	// in order to autonomously test/copmile this C++ file:
	#define DIMENSIONS 3
	namespace GRMHD { constexpr int nVar = 19; }
	int main() { return 0; }
#else
	// If you include to ExaHyPE instead:
	#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
	#include "AbstractGRMHDSolver_FV.h" // Defines:
	namespace GRMHD { constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables; }
#endif

#include "PDE.h"
#include "really.h"

void GRMHD::PDE::prepare() {
	// NEED to prepare:
	// Generic for Sij and preparation:
	DFOR(i) CONTRACT(j) Bmag.lo(i) = gam.lo(i,j) * Bmag.up(j);
	DFOR(i) CONTRACT(j) vel.lo(i)  = gam.lo(i,j) * vel.up(j);
	// S^i is needed in both flux and ncp
	DFOR(i) CONTRACT(j) Si.up(i)   = gam.up(i,j) * Si.lo(j);
	
	WW = SQ(Dens/rho); // W^2
	BmagBmag = 0; CONTRACT(k) BmagBmag += Bmag.lo(k)*Bmag.up(k); // B^2
	BmagVel = 0;  CONTRACT(j) BmagVel  += Bmag.up(j)*vel.lo(j); // B^j * v_j
	ptot = press + 0.5*(BmagBmag/WW + SQ(BmagVel)); // total pressure incl magn. field
}

void GRMHD::PDE::flux(double** Fluxes) {
	// Fluxes Shape: F[nDim][nVar]
	
	// Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
	Ul<sym::stored<3>> Sij;
	DFOR(i) DFOR(j) Sij.ul(i,j) = Si.up(i)*vel.lo(j) + ptot*delta(i,j) - Bmag.up(i)*Bmag.lo(j)/WW - BmagVel * vel.up(i) * Bmag.lo(j);
	
	// Zeta is the transport velocity (curly V in BHAC paper)
	Up<vec::stored<3>> zeta; DFOR(k) zeta.up(k) = alpha*vel.up(k) - beta.up(k);

	DFOR(k) { // F^k: flux in direction k.
		ConservedVariableShadow Flux(Fluxes[k]);
		Flux.Dens = Dens * zeta.up(k);
		DFOR(i) Flux.Si.lo(i) = alpha*Sij.ul(k,i) - beta.up(k)*Si.lo(i);
		Flux.tau = alpha*(Si.up(k) - vel.up(k)*Dens) - beta.up(k)*tau;
		DFOR(j) Flux.Bmag.up(j) = zeta.up(k)*Bmag.up(j) - zeta.up(j)*Bmag.up(k);
		
		// Constraint damping contributions:
		DFOR(j) Flux.Bmag.up(j) -= Bmag.up(k)*beta.up(j);
		Flux.phi = alpha*Bmag.up(k) - phi*beta.up(k);
	}
	
	// TODO :Zeroize the fluxes in the material parameters!
}

void GRMHD::PDE::RightHandSide(const double* const gradQ_Data, double* fusedSource_Data) {
	// This is the fusedSource = -NCP + AlgebraicSource

	ConservedVariableShadow Source(fusedSource_Data);
	const Gradients grad(gradQ_Data, nVar);
	
	// Sij is the 3-Energy-Momentum tensor. We need S^{ij} and S^i_j in the NCP.	
	UpSym<sym::stored<3>, sym::stored<3>> Sij;
	DFOR(i) DFOR(j) Sij.up(i,j) = Si.up(i)*vel.up(j) + gam.up(i,j)*ptot - Bmag.up(i)*Bmag.up(j)/WW - BmagVel * vel.up(i) * Bmag.up(j);
	Sij.ul_from_up(gam);

	// Source for D
	Source.Dens = 0;
	
	// U is the total energy seen by Eulerian observer
	// D is the density
	// tau is the rescaled energy density.
	// We evolve tau for an improved low density region.
	double U = tau + Dens;

	// S_i
	DFOR(i) {
		Source.Si.lo(i) = - U * grad.lo(i).alpha;
		CONTRACT(k) Source.Si.lo(i) += Si.lo(k) * grad.lo(i).beta.up(k);
		CONTRACT2(l,m) Source.Si.lo(i) += alpha/2. * Sij.up(l,m) * grad.lo(i).gam.lo(l,m);
	}
	
	// tau
	CONTRACT(k) Source.tau = -Si.up(k) * grad.lo(k).alpha; // Alejandros term
	// further BHAC terms for tau:
	CONTRACT2(i,j) Source.tau += Sij.ul(j,i) * grad.lo(j).beta.up(i);
	CONTRACT3(i,k,j) Source.tau += 1./2. * Sij.up(i,k)*beta.up(j) * grad.lo(j).gam.lo(i,k);
	
	// Starting from here:
	// Divergence Cleaning/Constraint damping sources for B^j and Phi.
	// In case no cleaning would be applied, Source.Bmag == Source.phi == 0.
	
	// B^i magnetic field
	DFOR(i) {
		Source.Bmag.up(i) = 0;
		CONTRACT(k) Source.Bmag.up(i) -= Bmag.up(k) * grad.lo(k).beta.up(i);
		CONTRACT(j) Source.Bmag.up(i) -= alpha * gam.up(i,j) * grad.lo(j).phi;
	}
	
	// phi
	Source.phi = - alpha * damping_term_kappa * phi; // algebraic source
	CONTRACT(k) Source.phi += Bmag.up(k) * grad.lo(k).alpha;
	CONTRACT(k) Source.phi -= phi * grad.lo(k).beta.up(k);
	CONTRACT3(k,l,m) Source.phi += alpha/2 * gam.up(l,m) * beta.up(k) * grad.lo(k).gam.lo(l,m);
}

/*
void GRMHD::PDE::algebraicSource(double* Source_data) {
	ConservedVariables Source(Source_data);
	// we could put the following term here:
	// Source.phi = - alpha * kappa * phi; // algebraic source
	// All other terms should be zero.
}

void GRMHD::PDE::fusedSource(const double* const gradQ_Data, double* Source) {
	// Variables Source(Source_data, nVar);
	// 1. Compute NCP
	nonConservativeProduct(gradQ_Data, Source);
	// 2. Flip sign
	NVARS(m) Source[m] = - Source[m];
	// 3. Add the algebraic source terms:
	// CONTRACT2(l,m) Source.tau += lapse * Sij.up(l,m) * Kextr.lo(l,m);
	// Source.phi = - lapse * kappa * phi;
	// 4. Done.
}
*/

void GRMHD::PDE::eigenvalues(const double* const Q, const int d, double* lambda) {
	//std::fill_n(lambda, nVar, 1.0);
	NVARS(m) lambda[m] = 1.0;
}
