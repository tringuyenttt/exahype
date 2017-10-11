/**
 * This is the plotter abused to compute global spacetime integrals.
 * This is actually used to compare with the exact solution for
 * convergence tests.
 **/

#include "PDE/PDE.h"
#include "Writers/IntegralsWriter.h"
#include "InitialData.h"
#include "kernels/GaussLegendreQuadrature.h"
#include <cmath>

#include "kernels/aderdg/generic/c/sizes.cpph"



GRMHD::IntegralsWriter::IntegralsWriter(exahype::solvers::LimitingADERDGSolver&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

GRMHD::IntegralsWriter::IntegralsWriter(GRMHD::GRMHDSolver_ADERDG&  solver)
	: IntegralsWriter() { plotForADERSolver = true; }

GRMHD::IntegralsWriter::IntegralsWriter(GRMHD::GRMHDSolver_FV&  solver)
	: IntegralsWriter() { plotForADERSolver = false; }

GRMHD::IntegralsWriter::IntegralsWriter() :
	conserved("output/cons-"),
	primitives("output/prim-"),
	errors("output/error-"),
	statistics("output/volform")
{
	conserved.add(0, "dens");
	conserved.add(1, "sconx");
	conserved.add(2, "scony");
	conserved.add(3, "sconz");
	conserved.add(4, "ener");
	conserved.add(5, "bx");
	conserved.add(6, "by");
	conserved.add(7, "bz");
	conserved.add(8, "psi");
		
	primitives.add(0, "rho"); // V[0]=Q[0]
	primitives.add(1, "velx");
	primitives.add(2, "vely");
	primitives.add(3, "velz");
	primitives.add(4, "press");
	primitives.add(5, "bx");
	primitives.add(6, "by");
	primitives.add(7, "bz");
	primitives.add(8, "psi");
	
	errors.add(0, "rho");
	errors.add(1, "velx");
	errors.add(2, "vely");
	errors.add(3, "velz");
	errors.add(4, "press");
	errors.add(5, "bx");
	errors.add(6, "by");
	errors.add(7, "bz");
	errors.add(8, "psi");
}


GRMHD::IntegralsWriter::~IntegralsWriter() {
	// delete all the TimeSeriesReductions. Not that important
	// as this is a kind of a Singleton object.
}


void GRMHD::IntegralsWriter::startPlotting(double time) {
	conserved.startRow(time);
	primitives.startRow(time);
	errors.startRow(time);
	statistics.startRow(time);
}


void GRMHD::IntegralsWriter::finishPlotting() {
	conserved.finishRow();
	primitives.finishRow();
	errors.finishRow();
	statistics.finishRow();
}


void GRMHD::IntegralsWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
	// make sure this plotter has no output associated
	assertion( outputQuantities == nullptr );


	double dV;
	if(plotForADERSolver) {
		const int order = GRMHD::AbstractGRMHDSolver_ADERDG::Order;
		dV = kernels::ADERDGVolume(order, sizeOfPatch, pos);
	} else {
		const int patchSize = GRMHD::AbstractGRMHDSolver_FV::PatchSize;
		dV = tarch::la::volume(sizeOfPatch)/patchSize; // correct is probably (patchSize+1)
	}

	statistics.addValue(dV, 1);


	// reduce the conserved quantities
	conserved.addValue(Q, dV);

	// reduce the primitive quantities
	double V[nVar];
	GRMHD::Cons2Prim(V, Q).copyFullStateVector();
	primitives.addValue(V, dV);

	// now do the convergence test, as we have exact initial data
	double ExactCons[nVar], ExactPrim[nVar];
	
	//id->Interpolate(xpos, timeStamp, ExactCons);
	AlfenWaveCons(x.data(),timeStamp,ExactCons);
	GRMHD::Cons2Prim(ExactPrim, ExactCons).copyFullStateVector();
	
	double localError[nVar];
	for(int i=0; i<nVar; i++) {
		localError[i] = std::abs(V[i] - ExactPrim[i]);
	}
	
	errors.addValue(localError, dV);
}


