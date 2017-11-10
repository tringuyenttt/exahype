
#include "InitialData/InitialData.h"
#include "tarch/logging/Log.h"


InitialDataCode& InitialDataCode::getInstance() {
	static tarch::logging::Log _log("InitialDataCode");
	static InitialDataCode* id = nullptr;
	if(id) return *id;

	if(1) {
		logInfo("InitialDataCode", "Loading AlfenWaveCons Initial Data");
		id = new AnalyticalID<AlfenWaveCons>;
	}

	if(0) {
		logInfo("getInstance()", "Loading PizzaTOV Initial Data");
		id = new pizzatov();
	}
	
	if(0) {
		logInfo("getInstance()", "Loading RNSID Initial Data");
		id = new rnsid();
	}
	
	return *id;
	// then, some day do the job of GlobalID.cpp:
	/*
	if(idname == "Fortran") {
		id = new fortranid();
		return true;
	} else if(idname == "PizzaTOV") {
		id = new pizzatov();
		return true;
	} else if(idname == "RNSID") {
		id = new rnsid();
		return true;
	}
	return false; // no success
	*/
}

#include "PDE/PDE.h"
#include "GRMHDSolver_ADERDG_Variables.h"
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

/// a shorthand
void InitialData(const double* x, double t, double* Q) {
	for(int i=0;i<nVar;i++) Q[i] = 0.0; // Zeroize
	
	InitialDataCode::getInstance().Interpolate(x,t,Q);
	
	  // also store the positions for debugging
	GRMHD::AbstractGRMHDSolver_ADERDG::Variables var(Q);
	GRMHD::AbstractGRMHDSolver_ADERDG::VariableShortcuts Qi;
	constexpr int posI = 19; // position 19
	for(int i=0; i<DIMENSIONS; i++) Q[posI+i] = x[i];
	if(DIMENSIONS<3) Q[posI+2] = -1;
	var.check() = 47110815;
	//zero2Din3D(Q);

	for(int i=0;i<nVar;i++) { if(!std::isfinite(Q[i])) { printf("Qid[%d] = %e\n", i, Q[i]); std::abort(); } }
}
