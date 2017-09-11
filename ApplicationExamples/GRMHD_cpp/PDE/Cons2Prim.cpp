#include "PDE.h"
#include <cmath> // max, sqrt
#include <algorithm> // max

#include <stdio.h> // printf

using namespace GRMHD;
using namespace std; // sqrt, max

// The C2P-Routines ported from Fortran for the SRMHD by Dominic, fully correct also for GRMHD.
double RTSAFE_C2P_RMHD1(const double X1, const double X2, const double XACC, const double gam, const double d, const double e, const double s2, const double b2, const double sb2, double& w,bool& failed);
void FUNC_C2P_RMHD1(const double x,double* f,double* df,const double gam,const double d,const double e,const double s2,const double b2, const double sb2,double& w);

// Remember: The removal or adding of \sqrt{gamma} to the quantities.

// At this stage, we can prepare conservative variables but do not yet have access to primitive ones.
void GRMHD::Cons2Prim::prepare() {
	// B_i: Needed for Sij and preparation:
	// v^i: is needed for the PDE system (Sij) and is either computed in the followup() or by the Cons2Prim operation.
	// S^i: is needed in both flux and ncp
	Bmag.lo=0; DFOR(i) CONTRACT(j) Bmag.lo(i) += gam.lo(i,j) * Bmag.up(j);
	Si.up  =0; DFOR(i) CONTRACT(j) Si.up(i)   += gam.up(i,j) * Si.lo(j);
	BmagBmag = 0; CONTRACT(k) BmagBmag += Bmag.lo(k)*Bmag.up(k); // B^j * B_j // needed for ptot
	SconScon = 0; CONTRACT(k) SconScon +=   Si.lo(k)*Si.up(k);   // S^j * S_j // needed for c2p
	BmagScon = 0; CONTRACT(k) BmagScon += Bmag.lo(k)*Si.up(k);   // B^j * S_j // needed for c2p
}

// At this stage, we have all primitives recovered and can postcompute some quantities. This is as a
// service or can be used for 
void GRMHD::Cons2Prim::followup() {
	BmagVel = 0;  CONTRACT(j) BmagVel  += Bmag.up(j)*vel.lo(j);  // B^j * v_j // needed for ptot
	WW = SQ(WLorentz); // W^2: Lorentz factor squared
	ptot = press + 0.5*(BmagBmag/WW + SQ(BmagVel)); // total pressure incl magn. field, needed in 3-energy-mom-tensor
}

void GRMHD::Cons2Prim::copyFullStateVector() {
	// 1) Assume that hydro variables have been set
	// 2) Copy only:
	copy_magneto(V);
	copy_adm(V);
}

void GRMHD::Cons2Prim::perform() {
	constexpr double tol       = 1e-8;
	constexpr double p_floor   = 1.0e-5;
	constexpr double rho_floor = 1.0e-4;
	
	constexpr double gamma = GRMHD::Parameters::gamma;
	
	// TODO here: Removal of 1./\sqrt{gamma}.

	// First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
	constexpr double gamma1    = (gamma - 1.0)/gamma;
	double e      = tau; // sic! My naming is probably sick.

	constexpr double eps    = 1.e-10;
	constexpr double x1     = 0.;
	constexpr double x2     = 1.-eps;

	// RTSAFE_C2P_RMHD1 has output {Gamma Factor w, Squared 3-velocity v^2}.
	bool failed=false;
	WLorentz=0;
	VelVel = RTSAFE_C2P_RMHD1(x1,x2,tol,gamma1,Dens,e,SconScon,BmagBmag,BmagScon*BmagScon,WLorentz,failed);
	
	if (failed) {
		// We should raise an error instead, the c2p failed.
		printf("C C2P FAILED\n");
		rho = rho_floor;
		press = p_floor;
		DFOR(i) vel.up(i) = 0;
	} else {
		double den  = 1.0/(WLorentz+BmagBmag);
		double vb   = BmagScon/WLorentz;
		rho  = Dens*sqrt(1.-VelVel);
		DFOR(i) vel.up(i) = (Si.up(i) + vb * Bmag.up(i))*den; // TODO: This looks wrong. CHECK
		DFOR(i) vel.lo(i) = (Si.lo(i) + vb * Bmag.lo(i))*den;
		press     = gamma1*(WLorentz*(1.-VelVel)-rho); // EOS
		press     = max(1.e-15, press); // bracketing
	}
}

double RTSAFE_C2P_RMHD1(const double X1, const double X2, const double XACC, const double gam, const double d,
    const double e, const double s2, const double b2, const double sb2, double& w,bool& failed) {
  int MAXIT=200;
  failed = false;
  double FL;
  double DF;

#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&X1,&FL,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(X1,&FL,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  double FL_test,DF_test,w_test;
  FUNC_C2P_RMHD1(X1,&FL_test,&DF_test,gam,d,e,s2,b2,sb2,w_test);
//  assertionNumericalEquals(FL,FL_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);

  double v2 = 0; // don't know if this is a good init value but flag failed should safe us
  if(FL==0) {
    v2=X1;
    return v2;
  }

  double FH;
#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&X2,&FH,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(X2,&FH,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  double FH_test;
  FUNC_C2P_RMHD1(X2,&FH_test,&DF_test,gam,d,e,s2,b2,sb2,w_test);
//  assertionNumericalEquals(FH,FH_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);

  if(FH==0) {
     v2=X2;
     return v2;
  }
  if(FL*FH>0) {
     failed = true;
     return v2;
  }

  double XL,XH;
  double SWAP;
  if(FL<0) {
    XL=X1;
    XH=X2;
  } else {
    XH=X1;
    XL=X2;

    SWAP=FL;
    FL=FH;
    FH=SWAP;
  }
  v2=.5*(X1+X2);
  double DXOLD=std::abs(X2-X1);
  double DX=DXOLD;

  double F;
#ifdef USE_FORTRAN_HELPER_FUNCS
  func_c2p_rmhd1_(&v2,&F,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
  FUNC_C2P_RMHD1(v2,&DF,&DF,gam,d,e,s2,b2,sb2,w);
#endif
  //double F_test;
  //FUNC_C2P_RMHD1(v2,&F_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
//  assertionNumericalEquals(F,F_test);
//  assertionNumericalEquals(DF,DF_test);
//  assertionNumericalEquals(*w,w_test);
  for (int J=1; J<MAXIT; ++J) {
     if(((v2-XH)*DF-F)*((v2-XL)*DF-F)>=0
          || std::abs(2.*F)>std::abs(DXOLD*DF) ) {
        DXOLD=DX;
        DX=0.5*(XH-XL);
        v2=XL+DX;
        if(XL==v2) {
          return v2;
        }
     } else {
        DXOLD=DX;
        DX=F/DF;
        double TEMP=v2;
        v2=v2-DX;
        if (TEMP==v2) {
          return v2;
        }
     }
     if (std::abs(DX)<XACC) {
       return v2;
     }
#ifdef USE_FORTRAN_HELPER_FUNCS
     func_c2p_rmhd1_(&v2,&F,&DF,&gam,&d,&e,&s2,&b2,&sb2,w);
#else
     FUNC_C2P_RMHD1(v2,&F,&DF,gam,d,e,s2,b2,sb2,w);
     //FUNC_C2P_RMHD1(v2,&F_test,&DF_test,gam,d,e,s2,b2,sb2,&w_test);
#endif
//     assertionNumericalEquals(F,F_test);
//     assertionNumericalEquals(DF,DF_test);
//     assertionNumericalEquals(*w,w_test);
     if(F<0) {
        XL=v2;
        FL=F;
     } else {
        XH=v2;
        FH=F;
     }
  }
  failed = true;
  return v2;
}


void FUNC_C2P_RMHD1(const double x,double* f,double* df,const double gam,const double d,const double e,const double s2,const double b2,
    const double sb2,double& w) {
  //
  // This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  // and it corresponds to their choice 3 in Section 3.2
  //
  double third = 1./3.;

  double v2=x;
  double rho=d*std::sqrt(1.-v2);

  double c3=1.-gam*(1.-v2);
  double c2=gam*rho+.5*b2*(1.+v2)-e;
  double c0=-0.5*sb2;

  //  For every x=v2, we solve a cubic in W of the form:
  //  c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  //  W=y of the paper. If sb=0 ( meaning c0 = 0),
  //  w = -c2/c3 > 0 and dw = 0 in the do loop below.
  //  If -c2/c3 < 0 when sb=0, which makes w=0,
  //  this is a signature that something was wrong before.
  if ( abs ( c0) < 1.0e-20) {
    w = -c2 / c3;
  } else {
    w = std::max ( - c2 / c3, std::pow(( -c0 / c3), third));   // ( -c0 / c3)**third
    for(int iter = 0; iter < 200; iter++) {
      double dw = -((c3*w + c2)*w*w + c0)/((3*c3*w + 2*c2)*w);
      if (std::abs(dw/w)<1.e-10) {
        break; // happy breakdown
      }
      w = w + dw;
    }
  };

  double dc3   = gam;
  double dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2));
  double dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2);
  double wb    = w + b2;
  double vb2 = sb2 / (w*w);
  *f   = wb*wb * v2 - ( 2.0 * w + b2) * vb2 - s2;
  *df  = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2));

//  std::cout << "w_out="<<*w_out;
//  std::cout << ",f="<<*f;
//  std::cout << ",df="<<*df << std::endl;
}
