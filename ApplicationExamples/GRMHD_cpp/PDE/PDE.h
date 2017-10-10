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
	//int main() { return 0; }
#else
	// If you include to ExaHyPE instead:
	#include "peano/utils/Dimensions.h" // Defines DIMENSIONS
	#include "AbstractGRMHDSolver_FV.h" // Defines: nVar
	namespace GRMHD { constexpr int nVar = GRMHD::AbstractGRMHDSolver_FV::NumberOfVariables; } // ensure this is 19 or so
	//namespace GRMHD { constexpr int nVar = 19; }
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

#include <cstring> // memcpy

//namespace GRMHD { constexpr int nVar = 42; }
#define NVARS(i) for(int i=0; i < GRMHD::nVar; i++)

namespace GRMHD {
	using namespace tensish; 
	using sym::delta; // Kronecker symbol
	
	namespace Parameters {
		// a temporary space for global system parameters
		
		// Ideal EOS:
		// 4/3 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
		// 2.0 used for TOV stars
		constexpr double gamma = 2.0;
      
		// Divergence cleaning:
		// 1.0 used in ADERDG3D-SRMHD-Z4 by MD, 01.Nov2016
		constexpr double DivCleaning_a = 1.0;
	}
	
	namespace Hydro {
		// Hydro vars are always with offset = 0.
		constexpr int size = 5;
		
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
				// Conserved hydro variables
				constexpr int Dens = 0;
				constexpr int Si_lo = 1;
				constexpr int tau = 4;
			}
			
			template<typename state_vector, typename vector_lo, typename scalar>
			struct StateVector {
				state_vector Q;
				scalar &Dens, &tau;
				vector_lo Si; /* sic */
				
				constexpr StateVector(state_vector Q_) :
					Q(Q_),
					Dens (Q[Indices::Dens]), // Conserved scalars
					tau  (Q[Indices::tau]),
					Si   (Q+Indices::Si_lo) // Conserved vectors
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
				const ConstLo<vec::const_shadow_D>,
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
				Lo<vec::shadow_D>,
				double> Shadow;
				
			/**
			* Access to the full set of conserved variables in a read only fashion which
			* however allows you to retrieve B_i and S^j, i.e. we have storage for these
			* variables. These two vectors are the only non-constant members of the
			* data structure as you need to fill them after construction.
			**/
			typedef StateVector<
				const double* const,
				UpLo<vec::stored_D, vec::const_shadow_D>::ConstLo,
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
				UpLo<vec::stored_D, vec::shadow_D>::InitLo,
				double> ShadowExtendable;

			// to be sorted:
			void toPrimitives(double* V); // C2P Cons2Prim operation
		} // ns ConservedHydro
		
		/**
		* Prmitive Variables.
		* They should only store the hydro variables, not Maxwell or ADM variables.
		**/
		namespace Primitives {
			namespace Indices {
				// The positions of the primitive hydro variables where there exists
				// a 1:1 mapping to the conserved hydro variables.
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
			};
		
			typedef StateVector<double*, Up<vec::shadow_D>, double> Shadow;
			typedef StateVector<const double* const, UpLo<vec::const_shadow_D,vec::stored_D>::ConstUp, const double> ConstShadowExtendable;
			typedef StateVector<double*, UpLo<vec::shadow_D, vec::stored_D>::InitUp, double> ShadowExtendable;
			struct Stored : public Shadow { double V[size]; Stored() : Shadow(V) {} };
		} // ns Primitives
	} // ns Hydro

	namespace Magneto {
		constexpr int size = 4; // length of Magneto/Maxwell-vector
		
		namespace RelativeIndices {
			constexpr int Bmag_up = 0;
			constexpr int Phi = 4;
		}
		namespace AbsoluteIndices {
			// Conserved (C2P invariant) Magnetic/Maxwell variables
			constexpr int offset = 5;
			constexpr int Bmag_up = 5;    /* length: 3 */
			constexpr int Phi = 8;
		}
		
		template<typename state_vector, typename vector_up, typename scalar>
		struct StateVector {
			state_vector Q_Magneto;
			scalar &phi;
			vector_up Bmag; /* sic */
			
			constexpr StateVector(state_vector Q) :
				Q_Magneto(Q),
				phi  (Q[AbsoluteIndices::Phi]),   // Conserved scalars
				Bmag (Q+AbsoluteIndices::Bmag_up) // Conserved vectors
				{}
			
			void copy_magneto(double* target) {
				std::memcpy(target + AbsoluteIndices::offset, Q_Magneto + AbsoluteIndices::offset, size*sizeof(double));
			}
		};
		
		typedef StateVector<const double* const, const ConstUp<vec::const_shadow_D>, const double> ConstShadow;
		typedef StateVector<double* const, Up<vec::shadow_D>, double> Shadow;
		typedef StateVector<const double* const, UpLo<vec::const_shadow_D, vec::stored_D>::ConstUp, const double> ConstShadowExtendable;
		typedef StateVector<double* const, UpLo<vec::shadow_D, vec::stored_D>::InitUp, double> ShadowExtendable;
	} // ns Magneto

	namespace ADMBase {
		constexpr int size = 10;
		
		namespace RelativeIndices {
			constexpr int lapse = 0;
			constexpr int shift_lo = 1;  /* length: 3 */
			constexpr int gam_lo = 4;    /* length: 6 */
		}
		namespace AbsoluteIndices {
			// The positions of the material parameters in the GRMHD state vector. There is an offset
			// due to the size of the hydro and maxwell variables before.
			constexpr int offset = 9;
			constexpr int lapse = 9;
			constexpr int shift_lo = 10;  /* length: 3 */
			constexpr int gam_lo = 13;    /* length: 6 */
			//constexpr int detg = 19;    // not yet
		}
		
		template<typename state_vector, typename vector_up, typename metric_lo, typename scalar>
		struct StateVector { // Material Parameters
			state_vector Q_ADM;
			scalar &alpha;//, &detg; // Scalars: Lapse, Determinant of g_ij
			vector_up beta; // Shift vector: (Conserved) Material parameter vector
			metric_lo gam;  // 3-Metric: (Conserved) Material parameter tensor
			
			constexpr StateVector(state_vector Q) :
				Q_ADM(Q),
				alpha(Q[AbsoluteIndices::lapse]),
				//detg (Q[AbsoluteIndices::detg]),
				beta (Q+AbsoluteIndices::shift_lo),
				gam  (Q+AbsoluteIndices::gam_lo)
				{}
				
			void copy_adm(double* target) {
				std::memcpy(target + AbsoluteIndices::offset, Q_ADM + AbsoluteIndices::offset, size*sizeof(double));
			}
			
			// set all ADM components to zero. This is unphysical (not identity) for real metric but
			// practical if gradients of the metric are stored.
			void zero_adm() {
				alpha = beta.up = gam.lo = 0;
			}
		};
		
		// A version of the ADMVariables where beta_lo and the upper metric
		// are not recovered.
		typedef StateVector<
			const double* const,
			const ConstUp<vec::const_shadow_D>, // shift
			const ConstLo<sym::const_shadow_D>, // metric
			const double
			> ConstShadow;
			
		// A writable shadowed ADMBase. Writeable is only useful for setting the initial data.
		// You can only set the lower components of the metric.
		typedef StateVector<
			double* const,
			Up<vec::shadow_D>,
			Lo<sym::shadow_D>,
			double
			> Shadow;
		
		// This is what the PDE needs: A working metric (upper and lower) but only upper shift.
		// This is read only (especially inside metric3).
		typedef StateVector<
			const double* const,
			// ConstUp_Lo<vec::const_shadow, vec::stored_D>, // if you ever need beta_lo
			ConstUp<vec::const_shadow_D>, // if you never need beta_lo
			metric3, // could store only parts of the metric here, too.
			const double
			> Full;
		
	} // ns ADM

	namespace GRMHDSystem { // need a better name
		constexpr int size = Hydro::size + Magneto::size + ADMBase::size;
		
		// A fully writable and shadowed state vector for the conserved variables
		struct Shadow : public Hydro::Conserved::Shadow, public Magneto::Shadow, public ADMBase::Shadow {
			constexpr Shadow(double* const Q) : Hydro::Conserved::Shadow(Q), Magneto::Shadow(Q), ADMBase::Shadow(Q) {}
		};
		
		struct ConstShadow : public Hydro::Conserved::ConstShadow, public Magneto::ConstShadow, public ADMBase::ConstShadow {
			constexpr ConstShadow(const double* const Q) : Hydro::Conserved::ConstShadow(Q), Magneto::ConstShadow(Q), ADMBase::ConstShadow(Q) {}
		};
		
		// other terms for these objects:
		typedef Shadow StateVector;
		typedef ConstShadow Gradient;
	}

	// A type storing the gradients of the conserved vector in one direction.
	// Since it does not make sense to
	// retrieve the lower/upper components of gradients, we don't even allocate storage or
	// provide conversion strategies.

	// Gradients in each direction: partial_i * something
	struct Gradients {
		const GRMHDSystem::Gradient dir[DIMENSIONS];
		// Access an element in some direction
		const GRMHDSystem::Gradient& operator[](int i) const { return dir[i]; }
		// Access an element, indicating that it's partial_i, not partial^i
		const GRMHDSystem::Gradient& lo(int i) const { return dir[i]; }
		
		#define dbl const double* const // shorthand
		
		#if DIMENSIONS == 3
		Gradients(dbl gradQ) : dir{gradQ+0,gradQ+nVar,gradQ+2*nVar} {}
		Gradients(dbl Qx, dbl Qy, dbl Qz) : dir{Qx,Qy,Qz} {}
		#elif DIMENSIONS == 2
		Gradients(dbl gradQ) : dir{gradQ+0,gradQ+nVar} {}
		Gradients(dbl Qx, dbl Qy) : dir{Qx,Qy} {}
		Gradients(dbl Qx, dbl Qy, dbl Qz) : Gradients(Qx,Qy) {} // for convenience
		#endif
		
		// we need proper shadowing for this
		//Gradients(const double* const gradQ[DIMENSIONS]) {}
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
	struct Cons2Prim : public Hydro::Conserved::ConstShadowExtendable, public Hydro::Primitives::ShadowExtendable, public Magneto::ConstShadowExtendable, public ADMBase::Full {
		Cons2Prim(double* const V, const double*const Q_) :
			Hydro::Conserved::ConstShadowExtendable(Q_),
			Hydro::Primitives::ShadowExtendable(V),
			Magneto::ConstShadowExtendable(Q_),
			ADMBase::Full(Q_)
			{ prepare(); perform(); followup(); }

		// Quantities for computing the energy momentum tensor
		// as well as the C2P.
		double	WLorentz,	///< Lorentz factor (Gamma)
			WW,		///< Squared Lorentz factor (Gamma^2)
			RhoEnthWW,	///< rho*h*Gamma^2 as given by rtsafe
			SconScon,	///< S^2 = S_i*S^i: Squared conserved momentum
			BmagBmag,	///< B^2 = B_i*B^i: Squared magnetic field
			BmagVel,	///< B_i*v^i: Magn field times three velocity
			BmagScon,	///< B_i*S^i: Magn field times squared cons momentum
			VelVel,		///< V^2 = v_i*v^i: Squared three velocity
			ptot;		///< ptot = p_hydro + p_mag: Total pressure
		
		/// Prepares the Conserved B_i, S^i as well as the quantities neccessary to compute
		/// The energy momentum tensor.
		void prepare();
		
		/// Computes the primitive variables with knowledge of all conserved/adm/local class variables.
		void perform();
		
		/// A modern interface to the RTSAFE C2P root finding procedure
		bool rtsafe(double& x, double& y);
		
		/// Compute the magnetic pressure and thelike. 
		void followup();
		
		/**
		 * Copy the full state vector V, not only the hydro variables. You want this if you
		 * want a "traditional" Cons2Prim in the Trento sense. Just call it as
		 *    Cons2Prim(Q,V).copyFullStateVector()
		 **/
		void copyFullStateVector();
		
		struct Stored;
	};

	/// A version which does not write to a shadowed storage but a local one
	struct Cons2Prim::Stored : public Cons2Prim {
		double V[Hydro::size];
		Stored(const double* const Q_) : Cons2Prim(V, Q_) {}
	};
	

	struct PDE : public Cons2Prim::Stored {
		typedef GRMHDSystem::Shadow Source;
		typedef GRMHDSystem::Shadow NCP;
		typedef GRMHDSystem::Shadow Flux;
		
		// Constraint damping constant
		static constexpr double damping_term_kappa = Parameters::DivCleaning_a;
		
		PDE(const double* const Q) : Cons2Prim::Stored(Q) {}

		// fluxes moved to own class in order to be able to compute them also in only one direction.
		
		//void flux(Flux& flux, int direction); // not that easy again because in each direction we have shared Sij and zeta.
		//Flux& Flux(double* Fk);
		//Fluxes& Fluxes(double** F);
		//Fluxes& Fluxes(double* Fx, double* Fy, double* Fz);
		// void flux(/* const double* const Q, */ double** Fluxes); // -> Moved to struct Fluxes.
		
		/************************* High level classy access ***************************/
		
		/// The BgradQ
		PDE& nonConservativeProduct(const Gradients& grad, NCP& ncp);
		/// Sets the algebraic source
		void algebraicSource(Source& Source_data);
		/// Adds the algebraic source (assumes Source=0 or similiar)
		void addAlgebraicSource(Source& Source_data);
		/// Computes the fusedSource = algebraicSource - NCP. This is the classical RightHandSide.
		void fusedSource(const Gradients& grad, Source& source);
		
		/// Eigenvalues: Currently trivially 1. We provide static access in order
		/// to avoid unneccessary boilerplate work.
		static void eigenvalues(const double* const Q, const int d,double* lambda);
		
		/*************************** Double pointer access ****************************/
		
		PDE& nonConservativeProduct(const double* const gradQ, double* ncp)
			{ Gradients g(gradQ); NCP n(ncp); return nonConservativeProduct(g, n); }
		PDE& nonConservativeProduct(const double* const Qx, const double* const Qy, const double* const Qz, double* ncp)
			{ Gradients g(Qx,Qy,Qz); NCP n(ncp); return nonConservativeProduct(g,n); }
		void algebraicSource(double* Source_data)
			{ Source s(Source_data); algebraicSource(s); }
		void fusedSource(const double* const gradQ, double* source)
			{ Gradients g(gradQ); Source s(source); fusedSource(g, s); }
		void fusedSource(const double* const Qx, const double* const Qy, const double* const Qz, double* source)
			{ Gradients g(Qx,Qy,Qz); Source s(source); fusedSource(g, s); }
			
		void fusedSource(const double* const gradQ, Source& source)
			{ Gradients g(gradQ); fusedSource(g, source); }
		void fusedSource(const double* const Qx, const double* const Qy, const double* const Qz, Source& source)
			{ Gradients g(Qx,Qy,Qz); fusedSource(g, source); }
		

	};
	
	// Compute a single flux in some direction
	struct FluxBase : public PDE {
		Mixed<sym::stored_D> Sij; ///< Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
		Up<vec::stored_D> zeta;   ///< Zeta is the transport velocity (curly V in BHAC paper)
		typedef GRMHDSystem::Shadow Flux;
		FluxBase(const double* const Q) : PDE(Q) { prepare(); }
		void prepare();
		void compute(Flux& flux, int k);
	};
	
	// Compute all fluxes at once
	struct Fluxes : public FluxBase {
		Flux F[DIMENSIONS];
		#if DIMENSIONS == 3
		Fluxes(double** F, const double* const Q) : FluxBase(Q), F{F[0],F[1],F[2]} {computeAll();}
		Fluxes(double* Fx, double* Fy, double* Fz, const double* const Q) : FluxBase(Q), F{Fx,Fy,Fz} {computeAll();}
		#elif DIMENSIONS == 2
		Fluxes(double** F, const double* const Q) : FluxBase(Q), F{F[0],F[1]} {computeAll();}
		Fluxes(double* Fx, double* Fy, const double* const Q) : FluxBase(Q), F{Fx,Fy} {computeAll();}
		Fluxes(double* Fx, double* Fy, double* Fz, const double* const Q) : Fluxes(Fx,Fy,Q) {} // for convenience
		#endif
		void computeAll() { DFOR(i) compute(F[i], i); }
		void zeroMaterialFluxes();
	};

	/*
	/// This is the algebraic conserved flux
	struct Fluxes : public PDE {
		typedef GRMHDSystem::Shadow Flux;
		Flux F[DIMENSIONS];
		Mixed<sym::stored_D> Sij; ///< Sij is the 3-Energy-Momentum tensor: We only need S^i_j in the flux.
		Up<vec::stored_D> zeta;   ///< Zeta is the transport velocity (curly V in BHAC paper)
		
		// other constructors:
		Fluxes(double** _F, const double* const Q) : PDE(Q), F{ _F[0], _F[1], _F[2] } { prepare(); DFOR(i) compute(i); }
		Fluxes(double* Fx, double* Fy, double* Fz, const double* const Q) : PDE(Q), F{Fx,Fy,Fz} { prepare(); DFOR(i) compute(i); }

		void prepare(); ///< compute Sij, zeta.
		void compute(int k);
		void zeroMaterialFluxes();
	};
	*/
	
	
	/// A first attempt. This definition is useful for initial data.
	struct Prim2Cons : public Hydro::Conserved::ShadowExtendable, public Magneto::ConstShadowExtendable, public Hydro::Primitives::ConstShadowExtendable, public ADMBase::Full {
		double W, enth, BmagBmag, BmagVel, VelVel;
		
		/**
		 * Attention: We map {Hydro::Primitives,ADMBase,Magneto} <-> V
		 *            and only the Hydro::Conserved <-> Q.
		 * This is good if you set everything to the primitives.
		 * However, if you have a mixed V-Q setting (i.e. B in Q, etc)
		 * then you need a different mapping.
		 **/
		Prim2Cons(double*const Q_, const double* const V) :
			Hydro::Conserved::ShadowExtendable(Q_),
			Magneto::ConstShadowExtendable(V),
			Hydro::Primitives::ConstShadowExtendable(V),
			ADMBase::Full(V)
			{ prepare(); perform(); }
		
		void prepare();
		void perform();
		
		void copyFullStateVector(); ///< See Cons2Prim copyFullStateVector
	};

} // namespace GRMHD
#endif /* GRMHD_PDE_CPP */
