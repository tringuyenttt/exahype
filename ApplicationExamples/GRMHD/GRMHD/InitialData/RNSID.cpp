#include "AbstractGRMHDSolver_ADERDG.h"
#include "GRMHDSolver_ADERDG_Variables.h"
#include "InitialData/InitialData.h"
#include "Fortran/PDE.h"

#ifndef RNSID_AVAILABLE

#include <stdlib.h>
#include <stdio.h>

rnsid::rnsid() {
        printf("Cannot call RNSID as not compiled with -DRNSID_AVAILABLE");
        abort();
}

void rnsid::Interpolate(const double* x, double t, double* Q) {}

#else /* RNSID_AVAILABLE */

#include "rnsid/rnsid.h"

rnsid::rnsid() {
	id = new RNSID::rnsid();
	
	// A TOV star
	id->axes_ratio = 1.0;
	id->rnsid_rho_min = 1e-10;
	
	// mapping for quantity vector
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts var;
	id->adm_idx.gxx = var.gij + 0;
	id->adm_idx.gxy = var.gij + 1;
	id->adm_idx.gxz = var.gij + 2;
	id->adm_idx.gyy = var.gij + 3;
	id->adm_idx.gyz = var.gij + 4;
	id->adm_idx.gzz = var.gij + 5;
	
	id->adm_idx.alp    = var.lapse;
	id->adm_idx.shift1 = var.shift + 0;
	id->adm_idx.shift2 = var.shift + 1;
	id->adm_idx.shift3 = var.shift + 2;
	
	// Note these are the primitives? ...
	id->hydro_idx.rho   = var.rho;
	id->hydro_idx.velx  = var.vel + 0;
	id->hydro_idx.vely  = var.vel + 1;
	id->hydro_idx.velz  = var.vel + 2;
	id->hydro_idx.press = var.E;
	
	id->Run();
}
	
void rnsid::Interpolate(const double* x, double t, double* Q) {
	constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;
	double V[nVar] = {0.0};
	id->Interpolate(x, V);
	
	// treatment of the atmostphere PROBABLY not done by RNSID
	const double atmo_rho = 1e-13;
	const double atmo_press = 1e-7;
	if(V[0] < atmo_rho) {
		V[0] = atmo_rho;
	}
	if(V[4] < atmo_press) {
		V[4] = atmo_press;
	}
	
	pdeprim2cons_(Q, V);
}

#endif /* RNSID_AVAILABLE */
