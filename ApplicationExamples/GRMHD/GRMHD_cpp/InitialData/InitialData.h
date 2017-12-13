#ifndef __INITIAL_DATA_ADAPTER_GRMHD_cpp__
#define __INITIAL_DATA_ADAPTER_GRMHD_cpp__

#include "PDE/PDE-GRMHD-ExaHyPE.h"
#include "Parameters/mexa.h"
#include <cmath> 
#include <string>

/// An Interface to initial data creation.
struct InitialDataCode {
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
	
	virtual void readParameters(const mexa::mexafile& parameters) {} // by default do nothing.
	
	/// This is the actual initial data registration
	static InitialDataCode* getInstanceByName(std::string name);
};

// a shorthand, also for debugging
void InitialData(const double* x, double t, double* Q);

// global initial data
struct GlobalInitialData {
	InitialDataCode* id;
	std::string name; ///< set by setIdByName
	
	GlobalInitialData(InitialDataCode* _id, std::string _name) : id(_id), name(_name) {}
	
	// Returns true in case of success
	bool setIdByName(const std::string& name);
	
	/// Gives a singleton instance
	static GlobalInitialData& getInstance();
	
	/// Read the parameters for the id.
	void readParameters(const mexa::mexafile& parameters);
	
	// as a shorthand
	static InitialDataCode& getInitialDataCode() { return *(getInstance().id); }
};

#include "PizzaTOV.h"
#include "RNSID.h"
#include "AnalyticalID.h"



#endif /* __INITIAL_DATA_ADAPTER_GRMHD_cpp__ */
