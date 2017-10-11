#include "InitialData.h"
#include "MySWESolver_Variables.h"
#include "MySWESolver.h"

#include <cmath>

using namespace std;

///// 2D /////

#ifdef Dim2

/*
 * Constant water height with "exponential" hump in bathymetry.  
 */
void SWE::ExpBreakProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h()  = 4.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = exp(-pow((x[0]-5),2)-pow((x[1]-5),2));
  } else {
    vars.h()  = 4.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = 0.0;
  }
}

/*
 * Constant water height with discontinuous jump in bathymetry.  
 */
void SWE::DamBreakProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h()  = 4.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = 1.0;
  } else {
    vars.h()  = 4.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = 0.0;
  }
}

/*
 * Sea at rest (steady state).  
 */
void SWE::SeaAtRestProblem(const double* const x,double* Q) {
  MySWESolver::Variables vars(Q);

  if((x[0] -5) *(x[0] -5) + (x[1] -5) *(x[1] -5) < 2) {
    vars.h()  = 3.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = 1.0;
  } else {
    vars.h()  = 4.0;
    vars.hu() = 0.0;
    vars.hv() = 0.0;
    vars.b()  = 0.0;
  }
}


#endif

void SWE::initialData(const double* const x,double* Q) {
  DamBreakProblem(x,Q);
}
