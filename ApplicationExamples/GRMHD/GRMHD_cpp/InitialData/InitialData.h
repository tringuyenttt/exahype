#ifndef __INITIAL_DATA_ADAPTER_GRMHD_cpp__
#define __INITIAL_DATA_ADAPTER_GRMHD_cpp__

#include "PDE/PDE-GRMHD-ExaHyPE.h"
#include <cmath> 
#include <string>

/// An Interface to initial data creation.
struct InitialDataCode {
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
	static InitialDataCode* getInstanceByName(const std::string& name);
};

// a shorthand, also for debugging
void InitialData(const double* x, double t, double* Q);

// global initial data
struct GlobalInitialData {
	InitialDataCode* id;
	GlobalInitialData(InitialDataCode* _id) : id(_id) {}
	
	// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Gives a singleton instance
	static GlobalInitialData& getInstance();
	
	// as a shorthand
	static InitialDataCode& getInitialDataCode() { return *(getInstance().id); }
};

#include "PizzaTOV.h"
#include "RNSID.h"
#include "AnalyticalID.h"



#endif /* __INITIAL_DATA_ADAPTER_GRMHD_cpp__ */
