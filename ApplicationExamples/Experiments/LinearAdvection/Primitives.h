#ifndef __EULER_PRIMITIVES__
#define __EULER_PRIMITIVES__

#include "GeneratedConstants.h"

// Ideal equation of state's Gamma
static const double eos_gamma = 1.4; // cf. fluxes, eigenvalues, initial data

// helper: Square of a.
inline double SQ(double a) { return (a)*(a); }

// another pseudo vector helper
inline double SQ3(const double* const v, int start=1) {
	double res=0;
	for(int i=0;i<3; i++) res += SQ(v[i+start]);
	return res;
}

/** The reason for this being in a header is pure lazyness and triviality of the C2P for Euler */
/*
inline void cons2prim(double* V, const double* Q) {
	// Euler c2p
	double p = (eos_gamma-1)*(Q[4] - 0.5 * SQ3(Q) / Q[0] );

	V[0] = Q[0]; // fluid density
	V[1] = Q[1] / Q[0]; // fluid velocities
	V[2] = Q[2] / Q[0];
	V[3] = Q[3] / Q[0];
	V[4] = p; // fluid pressure
}

inline void prim2con(double* Q, const double* V) {
	Q[0] = V[0]; // fluid density
	Q[1] = V[0] * V[1]; // momentum
	Q[2] = V[0] * V[2];
	Q[3] = V[0] * V[3];
	// total energy = internal energy + kinetic energy
	Q[4] = V[4] / (eos_gamma-1) + 0.5*V[0] * SQ3(V);
}

*/
// skip the cons2prim for debugging
inline void cons2prim(double* V, const double* Q) {
	for(int i=0; i<MY_NUMBER_OF_VARIABLES; i++)
		V[i] = Q[i];
}

	
inline void prim2con(double* Q, const double* V) {
	cons2prim(Q,V);
}




#endif /* __EULER_PRIMITIVES__ */
