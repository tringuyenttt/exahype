/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * Read PDE.h for further documentation.
 * Written at 2017-08-06 by SvenK.
 **/

#include "PDE.h"
#include <cmath> // NAN

/******************************************************************************
 * Status of the PDE system:
 *   It is now equal to the schemes of BHAC and ECHO, i.e. it seems to be valid
 *   for the cowling approximation but not neccessarily for a dynamical
 *   spacetime.
 * 
 * Important TODO: Stripping of the determinant of the metric is missing.
 *****************************************************************************/

#define S(x) printf(#x " = %e\n", x);
#define SI(x) S(x(0));S(x(1));S(x(2));

void GRMHD::FluxBase::prepare() {
	// Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
	SYMFOR(i,j) Sij.ul(i,j) = Si.up(i)*vel.lo(j) + ptot*delta(i,j) - Bmag.up(i)*Bmag.lo(j)/WW - BmagVel* vel.up(i) * Bmag.lo(j);
	
	// Zeta is the transport velocity (curly V in BHAC paper)
	DFOR(k) zeta.up(k) = alpha*vel.up(k) - beta.up(k);
}

void GRMHD::FluxBase::compute(Flux& flux, int k) {
	// F^k: flux in direction k.

	flux.Dens = Dens * zeta.up(k);
	DFOR(i) flux.Si.lo(i) = alpha*Sij.ul(k,i) - beta.up(k)*Si.lo(i);
	flux.tau = alpha*(Si.up(k) - vel.up(k)*Dens) - beta.up(k)*tau;
	DFOR(j) flux.Bmag.up(j) = zeta.up(k)*Bmag.up(j) - zeta.up(j)*Bmag.up(k);
	
	// Constraint damping contributions:
	DFOR(j) flux.Bmag.up(j) -= Bmag.up(k)*beta.up(j);
	flux.phi = alpha*Bmag.up(k) - phi*beta.up(k);
}

void GRMHD::Fluxes::zeroMaterialFluxes() {
	DFOR(k) {
		F[k].alpha = 0;
		F[k].beta.up = 0;
		F[k].gam.lo = 0;
	}
}
GRMHD::PDE& GRMHD::PDE::nonConservativeProduct(const Gradients& grad, Source& ncp) {
	// This is the NCP = matrixB * gradQ
	
	// Sij is the 3-Energy-Momentum tensor. We need S^{ij} and S^i_j in the NCP.	
	UpSym<sym::stored<3>, sym::stored<3>> Sij;
	SYMFOR(i,j) Sij.up(i,j) = Si.up(i)*vel.up(j) + gam.up(i,j)*ptot - Bmag.up(i)*Bmag.up(j)/WW - BmagVel * vel.up(i) * Bmag.up(j);
	Sij.ul=0; SYMFOR(i,j) CONTRACT(k) Sij.ul(i,j) += Sij.up(i,k) * gam.lo(j,k); // Sij^i_j = Sij^{ik} gam_{jk}

	// ncp for D
	ncp.Dens = 0;
	
	// U is the total energy seen by Eulerian observer
	// D is the density
	// tau is the rescaled energy density.
	// We evolve tau for an improved low density region.
	double U = tau + Dens;

	// S_i
	DFOR(i) {
		ncp.Si.lo(i) = U * grad.lo(i).alpha;
		CONTRACT(k) ncp.Si.lo(i) -= Si.lo(k) * grad.lo(i).beta.up(k);
		CONTRACT2(l,m) ncp.Si.lo(i) -= alpha/2. * Sij.up(l,m) * grad.lo(i).gam.lo(l,m);
	}
	
	// tau
	CONTRACT(k) ncp.tau = Si.up(k) * grad.lo(k).alpha;
	// tau: Cowling approximation instead of extrinsic curvature
	CONTRACT2(i,j) ncp.tau -= Sij.ul(j,i) * grad.lo(j).beta.up(i);
	CONTRACT3(i,k,j) ncp.tau -= 0.5 * Sij.up(i,k)*beta.up(j) * grad.lo(j).gam.lo(i,k);
	
	//CONTRACT3(k,i,j) { printf("%d,%d,%d: ",k,i,j); S(grad.lo(k).gam.lo(i,j)); }
	//S(ncp.tau);
	
	// Starting from here:
	// Divergence Cleaning/Constraint damping sources for B^j and Phi.
	// In case no cleaning would be applied, ncp.Bmag == ncp.phi == 0.
	
	// B^i magnetic field
	DFOR(i) {
		ncp.Bmag.up(i) = 0;
		CONTRACT(k) ncp.Bmag.up(i) += Bmag.up(k) * grad.lo(k).beta.up(i);
		CONTRACT(j) ncp.Bmag.up(i) += alpha * gam.up(i,j) * grad.lo(j).phi;
	}
	
	// phi
	ncp.phi = 0; // remember the algebraic source: alpha * damping_term_kappa * phi
	CONTRACT(k) ncp.phi -= Bmag.up(k) * grad.lo(k).alpha;
	CONTRACT(k) ncp.phi += phi * grad.lo(k).beta.up(k);
	CONTRACT3(k,l,m) ncp.phi += alpha/2 * gam.up(l,m) * beta.up(k) * grad.lo(k).gam.lo(l,m);
	
	return *this; // chainable
}

void GRMHD::PDE::algebraicSource(Source& source) {
	source.Dens = source.tau = source.Si.lo = source.Bmag.up = source.phi = 0;
	addAlgebraicSource(source);
}

void GRMHD::PDE::addAlgebraicSource(Source &source) {
	source.phi -= alpha * damping_term_kappa * phi; // algebraic source
	// CONTRACT2(l,m) Source.tau += lapse * Sij.up(l,m) * Kextr.lo(l,m);
}

void GRMHD::PDE::fusedSource(const Gradients& grad, Source& source) {
	// This is the fusedSource = -NCP + AlgebraicSource

	// 1. Compute NCP
	nonConservativeProduct(grad, source);
	// 2. Flip sign (move from lhs to rhs)
	NVARS(m) source.Q[m] = -source.Q[m];
	// 3. Add the algebraic source
	addAlgebraicSource(source);
}

void GRMHD::PDE::eigenvalues(const double* const Q, const int d, double* lambda) {
	//std::fill_n(lambda, nVar, 1.0);
	NVARS(m) lambda[m] = 1.0;
}
