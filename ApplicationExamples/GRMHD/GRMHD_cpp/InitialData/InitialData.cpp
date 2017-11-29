
#include "InitialData/InitialData.h"
#include "tarch/logging/Log.h"


InitialDataCode* InitialDataCode::getInstanceByName(const std::string& idname) {
	if(idname == "AlfenWave") 	return new AnalyticalID<AlfenWaveCons>;
	if(idname == "PizzaTOV")	return new pizzatov();
	if(idname == "RNSID") 		return new rnsid();
	return nullptr;
}

GlobalInitialData& GlobalInitialData::getInstance() {
	static GlobalInitialData* me = nullptr;
	if(!me) me = new GlobalInitialData(nullptr);
	return *me;
}

bool GlobalInitialData::setIdByName(const std::string& name) {
	static tarch::logging::Log _log("InitialDataCode");
	if(id) logWarning("setIdByName()", "Have already set global initial data, now overwriting with "<< name << ".");
	id = InitialDataCode::getInstanceByName(name);
	if(id) {
		logInfo("setIdByName()", "Successfully loaded "<< name << " Initial Data.");
		return true;
	} else {
		logError("setIdByName()", "Requested Initial Data '"<< name << "' not known.");
		return false;
	}
}

#include "PDE/PDE.h"
#include "GRMHDSolver_ADERDG_Variables.h"
constexpr int nVar = GRMHD::AbstractGRMHDSolver_ADERDG::NumberOfVariables;

/// a shorthand
void InitialData(const double* x, double t, double* Q) {
	for(int i=0;i<nVar;i++) Q[i] = 0.0; // Zeroize
	
	if(GlobalInitialData::getInstance().id == nullptr) {
		static tarch::logging::Log _log("InitialData");
		logError("InitialData()", "Cannot access InitialData because no initial Data has been defined yet.");
		std::abort();
	}
	
	GlobalInitialData::getInitialDataCode().Interpolate(x,t,Q);
	
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
