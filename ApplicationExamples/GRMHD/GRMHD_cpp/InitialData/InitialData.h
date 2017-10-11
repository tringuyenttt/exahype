#ifndef __INITIAL_DATA_ADAPTER_GRMHD_cpp__
#define __INITIAL_DATA_ADAPTER_GRMHD_cpp__

#include "PDE/PDE.h"
#include <cmath> 
#include <string>

/// An Interface to initial data creation.
struct InitialDataCode {
	virtual void Interpolate(const double* x, double t, double* Q) = 0;
	
	/// get the global unique instance of an initial data solver.
	static InitialDataCode& getInstance();
};

// a shorthand, also for debugging
void InitialData(const double* x, double t, double* Q);

#include "PizzaTOV.h"
#include "RNSID.h"
#include "AnalyticalID.h"



#endif /* __INITIAL_DATA_ADAPTER_GRMHD_cpp__ */
