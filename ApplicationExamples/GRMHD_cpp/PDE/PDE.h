/**
 * A C++ variant of the GRMHD equations for ExaHyPE.
 * 
 * This formulation uses the Tensish library for the Tensor abstraction.
 * We mimic the ExaHyPE notation of state vector abstractions with structures.
 * We copy all neccessary parts for practical purposes (scalar, vector, tensor
 * notation) as well as distinguation between conserved, primitive, material
 * parameters.
 *
 * Written at 2017-08-06 by SvenK.
 **/
#ifndef GRMHD_PDE_CPP
#define GRMHD_PDE_CPP

#ifdef TEST_NEW_PDE_AUTONOMOUSLY
	// in order to autonomously test/copmile this C++ file:
	#define DIMENSIONS 3
	namespace GRMHD { constexpr int nVar = 19; }
	int main() { return 0; }
#else
	// If you include to ExaHyPE instead:
	#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
	#include "AbstractGRMHDSolver_FV.h" // Defines:
	namespace GRMHD { constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables; }
#endif

/* This header file needs the definition of a constexpr int GRMHD::nVar
 * as well as a preprocessor variable DIMENSIONS.
 * You should provide this in the including file, for instance:
 *
 *   namespace GRMHD { constexpr int nVar = 10; }
 *   #define DIMENSIONS 3
 *   #include "PDE.h"
 *
 * DIMENSIONS has to be a C preprocessor variable for C macro ifdefs.
 * As ExaHyPE, there are some places where this code only works for 2 and 3
 * dimensions.
 */

#include "tensish.cpph"

//namespace GRMHD { constexpr int nVar = 42; }
#define NVARS(i) for(int i=0; i < GRMHD::nVar; i++)

namespace GRMHD {
	using namespace tensish; 
	using sym::delta; // Kronecker symbol
	
	/**
	* Storage for the conserved MHD variables, i.e. the tuple
	*   Q = (Dens, S_i, tau, B^i, phi).
	* The structure allows addressing the vector Q with names.
	* The template allows to specify types for the scalars and 3-vectors inside Q.
	* See the typedefs below for usable instances of the "Conserved" structure.
	*
	* In principle, this structure mimics the ExaHyPE generated variables structure
	* but in a much more usable manner (exploiting all the substructure of Q).
	**/
	namespace Conserved {
		/**
		* These structs store the index positions of fields, i.e. Q.Dens=0.
		* They are similiar to the VariableShortcuts which are generated
		* by the ExaHyPE java toolkit. However, by defining our own we are
		* independent and can make use of the functional notation.
		**/
		namespace Indices {
			// Conservatives
			constexpr int Dens = 0;
			constexpr int Si_lo = 1;
			constexpr int tau = 4;
			// Conserved (C2P invariant)
			constexpr int Bmag_up = 5;    /* length: 3 */
			constexpr int Phi = 8;
		};
		
		template<typename state_vector, typename vector_lo, typename vector_up, typename scalar>
		struct StateVector {
			state_vector Q; // const double* const
			
			// Conserved Scalars
			scalar &Dens, &tau, &phi;
			
			// Conserved vectors
			vector_lo Si; /* sic */
			vector_up Bmag; /* sic */
			
			constexpr StateVector(state_vector Q_) :
				Q(Q_),
				// Conserved scalars
				Dens (Q[Indices::Dens]),
				tau  (Q[Indices::tau]),
				phi  (Q[Indices::Phi]),
				// Conserved vectors
				Si   (Q+Indices::Si_lo),
				Bmag (Q+Indices::Bmag_up)
				{}
		};
		
		/**
		* Access the Conserved Variables Q fully read-only without access to anything
		* beyond what's inside Q (no further storage allocated except pointers to Q).
		* That is, you cannot access S^i or B_i as we don't compute it. Therefore, we
		* say this structure is just "shadowing" Q. You might want to use it like
		* 
		*    const ConservedConstShadow Q(Q_);
		*    double foo = Q.Bmag.up(1) * Q.Dens;
		* 
		**/
		typedef StateVector<
			const double* const,
			const ConstLo<vec::const_shadow>,
			const ConstUp<vec::const_shadow>,
			const double> ConstShadow;
			
		/**
		* The variable version of the ConservedConstShadow: Allows to change all
		* entries of Q by their aliases, as well as directly changing the vectors, i.e.
		*
		*    ConservedVariableShadow Q(Q_);
		*    Q.Bmag.up(2) = Q.Dens + 15;
		* 
		* But not:
		* 
		*    Q.Bmag.lo(2) = 15; 
		* 
		* That is, we dont' allow setting (or even getting) B_i or S^j. This can be
		* helpful to avoid errors in the formulation.
		**/
		typedef StateVector<
			double* const,
			Lo<vec::shadow>,
			Up<vec::shadow>,
			double> Shadow;
			
		/**
		* Access to the full set of conserved variables in a read only fashion which
		* however allows you to retrieve B_i and S^j, i.e. we have storage for these
		* variables. These two vectors are the only non-constant members of the
		* data structure as you need to fill them after construction.
		**/
		typedef StateVector<
			const double* const,
			UpLo<vec::stored<3>,vec::const_shadow>::ConstLo,
			UpLo<vec::const_shadow, vec::stored<3>>::ConstUp,
			const double> ConstShadowExtendable;
			
		/**
		* The variable version of ConservedConstFull: You can change all conserved
		* variables as well access the B_i and S^j. Note that
		* 
		*   Q.Bmag.lo(2) = 15;
		* 
		* will not change the shadowed double* Q_.
		**/
		typedef StateVector<
			double* const,
			UpLo<vec::stored<3>, vec::shadow>::InitLo,
			UpLo<vec::shadow, vec::stored<3>>::InitUp,
			double> ShadowExtendable;

		// to be sorted:
		void toPrimitives(double* V); // C2P Cons2Prim operation
	} // ns ConservedHydro
	
	/**
	 * Prmitive Variables.
	 * They should only store the hydro variables.
	 **/
	namespace Primitives {
		namespace Indices {
			constexpr int rho = 0;
			constexpr int vec_up = 1; /* sic */
			constexpr int press = 4;
		}
		
		template<typename state_vector, /*typename vector_lo, */typename vector_up, typename scalar>
		struct StateVector {
			state_vector V;
			
			// Primitive Scalars
			scalar &rho, &press;
			
			// Primitive vectors
			vector_up vel;
			
			constexpr StateVector(state_vector V_) :
				V(V_),
				// Primitive Scalars
				rho  (V[Indices::rho]),
				press(V[Indices::press]),
				// Primitive vectors
				vel  (V+Indices::vec_up)
				{}
			
			void toConserved(double* Q); // P2C Prims2Cons operation
		};
	
		typedef StateVector<
			double*,
			Up<vec::shadow>,
			double
			> Shadow;
			
		typedef StateVector<
			const double* const,
			UpLo<vec::const_shadow, vec::stored<3>>::ConstUp,
			const double
			> ConstShadowExtendable;
		
		// Writable primitives with velocity up. dont know if we need that yet
		typedef StateVector<
			double*,
			UpLo<vec::shadow, vec::stored<3>>::InitUp,
			double
			> ShadowExtendable;
			
		struct Stored : public ShadowExtendable {
			double V[nVar];
			Stored() : ShadowExtendable(V) {}
		};

		// this is what PDE inherits from. It immediately does the cons2prim
		// => CANNOT WORK because we need the ADM variables and several reconstructions
		/*
		struct StoredFromConserved : public Stored {
			StoredFromConserved(const double* const Q_) : Stored() {
				// do the Cons2Prim here.
				// i.e. use Q -> write to V
			}
		};*/
	} // ns Primitives
		

	namespace ADMBase {
		namespace Indices {
			// material parameters
			constexpr int lapse = 9;
			constexpr int shift_lo = 10;  /* length: 3 */
			constexpr int gam_lo = 13;    /* length: 6 */
			constexpr int detg = 19;
		}
		
		template<typename state_vector, typename vector_up, typename metric_lo, typename scalar>
		struct StateVector { // Material Parameters
			scalar &alpha, &detg; // Scalars: Lapse, Determinant of g_ij
			vector_up beta; // Shift vector: (Conserved) Material parameter vector
			metric_lo gam;  // 3-Metric: (Conserved) Material parameter tensor
			
			constexpr StateVector(state_vector Q) :
				alpha(Q[Indices::lapse]),
				detg (Q[Indices::detg]),
				beta (Q+Indices::shift_lo),
				gam  (Q+Indices::gam_lo)
				{}
		};
		
		// A version of the ADMVariables where beta_lo and the upper metric
		// are not recovered.
		typedef StateVector<
			const double* const,
			const ConstUp<vec::const_shadow>, // shift
			const ConstLo<sym::const_shadow>, // metric
			const double
			> ConstShadow;
			
		// A writable shadowed ADMBase. Writeable is only useful for setting the initial data.
		// You can only set the lower components of the metric.
		typedef StateVector<
			double* const,
			Up<vec::shadow>,
			Lo<vec::shadow>,
			double
			> Shadow;
		
		// This is what the PDE needs: A working metric (upper and lower) but only upper shift
		typedef StateVector<
			const double* const,
			// ConstUp_Lo<vec::const_shadow, vec::stored<3>>, // if you ever need beta_lo
			ConstUp<vec::const_shadow>, // if you never need beta_lo
			metric3, // could store only parts of the metric here, too.
			const double
			> Full;
		
	} // ns ADM
	
	// A full state vector, containing the primitives and the ADMBase, as
	// read only shadowed option
	/// -> Has ambiguity whether we want to recover beta_lo, upper metric, etc.
	
	namespace HydroBase {
		// for convenience: A GRMHD HydroBase.
		
		// A fully writable and shadowed state vector for the conserved variables
		struct StateVector : public Conserved::Shadow, public ADMBase::Shadow {
			constexpr StateVector(double* const Q) : Conserved::Shadow(Q), ADMBase::Shadow(Q) {}
		};
		
		// etc.
	}

	// A type storing the gradients of the conserved vector in one direction.
	// Since it does not make sense to
	// retrieve the lower/upper components of gradients, we don't even allocate storage or
	// provide conversion strategies.
	struct Gradient : public ADMBase::ConstShadow, public Conserved::ConstShadow {
		constexpr Gradient(const double* const Q) : ADMBase::ConstShadow(Q), Conserved::ConstShadow(Q) {}
	};

	// Gradients in each direction: partial_i * something
	struct Gradients {
		const Gradient dir[DIMENSIONS];
		// Access an element in some direction
		const Gradient& operator[](int i) const { return dir[i]; }
		// Access an element, indicating that it's partial_i, not partial^i
		const Gradient& lo(int i) const { return dir[i]; }
		
		Gradients(const double* const gradQ, const int nVar) :
			#if DIMENSIONS == 3
			dir{gradQ+0, gradQ+nVar, gradQ+2*nVar}
			#elif DIMENSIONS == 2
			dir{gradQ+0, gradQ+nVar}
			#endif
			// function body:
			{}
	};

	/**
	 * The PDE structure wraps ("shadows") the read-only conserved vector Q in order to provide
	 * the possibility to compute the flux and the source terms. It inherits the conserved
	 * variables, primitive variables and ADM variables for convenient access when implementing
	 * the source and the flux.
	 *
	 * It used to provide storage for the 3-Energy momentum tensor in the past, but now this is
	 * stored directly in the methods because they need different forms (S_ij, S^ij or S^i_j)
	 * which we address differently in the current formalism.
	 **/
	struct PDE : public Conserved::ConstShadowExtendable, public Primitives::Stored, public ADMBase::Full {
		// 3-Energy momentum tensor
		// Full<sym::stored<3>, sym::stored<3>> Sij;
		
		// Constraint damping constant
		const double damping_term_kappa = 5;
		
		PDE(const double*const Q_) :
			Conserved::ConstShadowExtendable(Q_),
			Primitives::Stored(/*Q_*/),
			ADMBase::Full(Q_)
			{ prepare(); Cons2Prim(); }

		// Quantities for computing the energy momentum tensor
		double WW, SconScon, BmagBmag, BmagVel, ptot;
		
		/// Prepares the Conserved B_i, S^i as well as the quantities neccessary to compute
		/// The energy momentum tensor.
		void prepare();
		
		/// Computes the primitive variables with knowledge of all conserved/adm/local class variables.
		void Cons2Prim();
		
		/// This is the algebraic conserved flux. Chainable.
		void flux(/* const double* const Q, */ double** Fluxes);
		
		/// This is the fusedSource, but we just call it RightHandSide because it is on the RHS
		/// of the PDE. Chainable.
		void RightHandSide(/* const double* const Q, */const double* const gradQ_Data, double* Source_data);
		
		// You can recover these functions by just reworking what's in the fusedSource.
		//void nonConservativeProduct(/* const double* const Q, */ const double* const gradQ_Data, double* BgradQ_Data);
		//void algebraicSource(/* const double* const Q, */ double* Source_data);
		
		// currently trivial. Since we don't need to access Q, we provide static
		// access in order to avoid unneccessary boilerplate work.
		static void eigenvalues(const double* const Q, const int d,double* lambda);
	};

} // namespace GRMHD
#endif /* GRMHD_PDE_CPP */
