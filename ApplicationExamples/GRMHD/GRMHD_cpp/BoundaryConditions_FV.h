/* Outsourcing from GRMHDSolver_*.cpp */
#ifndef GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL
#define GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL

#include "GRMHDSolver_FV_Variables.h" // just for GRMHD-specific variable names
#include "InitialData/InitialData.h"
#include "PDE/PDE.h"
#include "Parameters/RuntimeParameters.h"


#define FV_BOUNDARY_SIGNATURE \
    const double* const x, const double t,const double dt,\
    const int faceIndex, const int d,\
    const double* const stateIn, double* stateOut
#define FV_BOUNDARY_CALL \
   x, t, dt, faceIndex, d, stateIn, stateOut

// The difference between the ADERDG and FV boundary signature is that in ADERDG, also
// the fluxIn and fluxOut pointers are given.
// Thus the ADERDG BC signature is a real superset of the FV BC signature and thus we
// stick to the ADERDG signature as the more generic one.

#define BOUNDARY_SIGNATURE   FV_BOUNDARY_SIGNATURE
#define BOUNDARY_CALL        FV_BOUNDARY_CALL


// Note that passing a number of variables like this is an ugly way to circumvent a proper
// object oriented structure which is passed instead.

template <typename SolverType>
class BoundaryConditions {
public:
	// face indices: 0 x=xmin 1 x=xmax, 2 y=ymin 3 y=ymax 4 z=zmin 5 z=zmax
	// corresponding 0-left, 1-right, 2-front, 3-back, 4-bottom, 5-top
	static constexpr int EXAHYPE_FACE_LEFT = 0;
	static constexpr int EXAHYPE_FACE_RIGHT = 1;
	static constexpr int EXAHYPE_FACE_FRONT = 2;
	static constexpr int EXAHYPE_FACE_BACK = 3;
	static constexpr int EXAHYPE_FACE_BOTTOM = 4;
	static constexpr int EXAHYPE_FACE_TOP = 5;
	
	SolverType* solver;
	
	// a callback to SolverType. Typically this is fulfilled by ADERDG::adjustPointSolution
	// or FV::adjustSolution (for inconsistency-reasons wrongly labeled).
	typedef void (SolverType::*adjustSolutionType)(double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,
		const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt);
	
	// constexpr SolverType::VariableShortcuts vars; // this does not work
	//const GRMHD::AbstractGRMHDSolver_FV::VariableShortcuts vars;
	
	typedef BoundaryConditions<SolverType> me;
	typedef void (me::*boundarymethod)(BOUNDARY_SIGNATURE);
	
	// this could be improved by using boundarymethod side[6]
	// and having apply defined just as  side[faceIndex](...);.
	boundarymethod left, right, front, back, bottom, top;
	
	BoundaryConditions(SolverType* _solver) : solver(_solver),
		left(nullptr), right(nullptr), front(nullptr), back(nullptr), bottom(nullptr), top(nullptr) {}
	
	bool allFacesDefined() {
		return (left != nullptr && right != nullptr && front != nullptr && back != nullptr && bottom != nullptr && top != nullptr);
	}
	
	/// Apply by looking up the saved boundarymethods
	void apply(BOUNDARY_SIGNATURE) {
		switch(faceIndex) {
			case EXAHYPE_FACE_LEFT:   (this->*left  )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_FRONT:  (this->*front )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BOTTOM: (this->*bottom)(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_RIGHT:  (this->*right )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_BACK:   (this->*back  )(BOUNDARY_CALL); break;
			case EXAHYPE_FACE_TOP:    (this->*top   )(BOUNDARY_CALL); break;
			default: throw std::runtime_error("Inconsistent face index");
		}
	}
	
	// just to highlight that apply() is the main method:
	// void operator()(ADERDG_BOUNDARY_SIGNATURE) { apply(ADERDG_BOUNDARY_CALL); }
	
	/***************************************************************************
	 * Boundary Methods.
	 * 
	 * To add your own, do the following:
	 * 
	 *  1) Implement it with a method void yourName(ADERDG_BOUNDARY_SIGNATURE)
	 *  2) Add yourName to the parseFromString thing
	 * 
	 ***************************************************************************/
	
	#define BC_FROM_SPECFILE_SIDE(name) \
		if(constants.is_string(#name)) { \
			std::string val = constants.get_string(#name);\
			name = parseFromString(val);\
			if(name==nullptr) {\
				logError("setFromSpecFile", "Boundary condition at " #name ": Invalid method: '" << val << "'");\
			} else {\
				logInfo("setFromSpecFile", #name " boundary: " << val);\
			}\
		} else {\
			logError("setFromSpecFile", "Boundary condition at " #name " is not valid according to specfile parser.");\
		}

	/**
	 * MapView shall have methods
	 *    bool isValueValidString(std::string)
	 *    std:.string getvalueAsString(std::string)
	 * such as the Parser::ParserView class has.
	 * 
	 * @returns true in case of success, false otherwise
	 **/
	bool setFromMapView(const RuntimeParameters::MapView& constants) {
		tarch::logging::Log _log("BoundaryCondition");
		BC_FROM_SPECFILE_SIDE(left);
		BC_FROM_SPECFILE_SIDE(right);
		BC_FROM_SPECFILE_SIDE(front);
		BC_FROM_SPECFILE_SIDE(back);
		BC_FROM_SPECFILE_SIDE(bottom);
		BC_FROM_SPECFILE_SIDE(top);
		return allFacesDefined();
	}

	/**
	 * assign a boundary method by some string value
	 * @returns pinter if value could be correctly parsed, in case of error nullptr.
	 **/
	static boundarymethod parseFromString(const std::string& value) {
		boundarymethod target=nullptr;
		if(value == "zero")
			target = &me::vacuum;
		if(value == "reflective" || value == "refl")
			target = &me::reflective;
		if(value == "copy" || value == "outflow")
			target = &me::outflow;
		if(value == "exact")
			target = &me::exact;
		return target;
		//return !(target==nullptr);
		//throw std::exception("Invalid boundary method"); // Todo: return a usable error, or std::optional<BoundaryMethod> or so.
	}

	static constexpr int nVar = SolverType::NumberOfVariables;
	static constexpr int nDim = DIMENSIONS;

	/**
	 * The vacuum is problem dependent, thus this is really MHD vacuum in
	 * Minkowski spacetime for GRMHD.
	 **/
	void vacuum(BOUNDARY_SIGNATURE) {
		VacuumInitialData foo(stateOut);
	}
	
	/// Outflow / Copy BC
	void outflow(BOUNDARY_SIGNATURE) {
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
		}
	}
	
	/// Reflection/Hard wall BC, Currently WRONG.
	void reflective(BOUNDARY_SIGNATURE) {
		// First, copy state
		for(int m=0; m<nVar; m++) {
			stateOut[m] = stateIn[m];
		}
		
		SVEC::GRMHD::Shadow Q(stateOut);
		
		// vectors: flip sign
		Q.Si.lo(d) *= -1;
		Q.Bmag.up(d) *= -1;
		
		// Note: ADM BC should never be used anywhere as
		// ADM variables are parameters in this setup anyway!
		
		Q.beta.up(d) *= -1;
		// Symmetric Tensors: flip sign of off-diagonal
		DFOR(i) if(i != d) Q.gam.lo(i,d) *= -1;
	}
		
	/// Exact BC
	void exact(BOUNDARY_SIGNATURE) {
		// FV: No time integration
		InitialData(x, t, stateOut);
	}
};

typedef BoundaryConditions<GRMHD::GRMHDSolver_FV> BoundaryConditionsFV;

#endif /* GENERIC_BOUNDARY_CONDITIONS_WHICH_SHOULD_GO_SOMEWHERE_CENTRAL */
