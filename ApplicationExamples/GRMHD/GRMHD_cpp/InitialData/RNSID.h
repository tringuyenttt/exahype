#ifndef EXAHYPE_GRMHD_RNSID
#define EXAHYPE_GRMHD_RNSID

#include "InitialData/InitialData.h"

// Forward declaration for the actual RNSID implementation.
namespace RNSID {
	class rnsid;
}// ns RNSID

/**
 * The RNSID initial data are only accessible if support is compiled into the
 * ExaHyPE binary.
 * 
 * This class uses the pimpl mechanism and has a stub in the RNSID.cpp implementation
 * in order to allow a seamless compilation in any way.
 **/
struct rnsid : public InitialDataCode {
	RNSID::rnsid *id;
	
	rnsid();
	void Interpolate(const double* x, double t, double* Q);
};

#endif /* EXAHYPE_GRMHD_RNSID */
