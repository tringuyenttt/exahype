#include "MyElasticSolver.h"

#include "MyElasticSolver_Variables.h"

void Elastic::MyElasticSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool Elastic::MyElasticSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) {
  // @todo Please implement/augment if required
   if (tarch::la::equals(t, 0.0, 1e-15)) {  // @todo precision
    return true;
  }
   return false;
}

void Elastic::MyElasticSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  // State variables:
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
}

exahype::solvers::Solver::RefinementControl Elastic::MyElasticSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Elastic::MyElasticSolver::eigenvalues(const double* const Q,const int normalNonZeroIndex,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
   double cs = 3.464;  // km/s
  double cp = 6.0;    // km/s
  
  lambda[0] = -cp;
  lambda[1] = -cs;
  lambda[2] = 0.0;
  lambda[3] = cs;
  lambda[4] = cp;
}

void Elastic::MyElasticSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;
}


void Elastic::MyElasticSolver::algebraicSource(const double* const Q,double* S) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}


void Elastic::MyElasticSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  // @todo Please implement/augment if required
  stateOut[0] = 0.0;
  stateOut[1] = 0.0;
  stateOut[2] = 0.0;
  stateOut[3] = 0.0;
  stateOut[4] = 0.0;

  fluxOut[0] = 0.0;
  fluxOut[1] = 0.0;
  fluxOut[2] = 0.0;
  fluxOut[3] = 0.0;
  fluxOut[4] = 0.0;
}


void Elastic::MyElasticSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  // @todo Please implement/augment if required

  double lam;
  double mu;
  double rho;

  double cs = 3.464;  // km/s
  double cp = 6.0;    // km/s

  rho = 2.67;   // gm/cm^3
  
  lam = rho * (cp * cp - 2. * cs * cs);  // lambda (MPa)
  mu= rho * cs * cs;

  // @todo Please implement/augment if required
  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  BgradQ[4] = 0.0;

  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = 0.0;
  BgradQ[8] = 0.0;
  BgradQ[9] = 0.0;

  const double* Qx = &gradQ[0]; // (:,1)
  const double* Qy = &gradQ[5]; // (:,2)
 


  BgradQ[0] = -1.0/rho*Qx[2];
  BgradQ[1] = -1.0/rho*Qx[4];

  BgradQ[2] = -(2.0*mu + lam)*Qx[0];
  BgradQ[3] = -lam*Qx[0];
 
  BgradQ[4] = -mu*Qx[1];
  

  BgradQ[5] = -1.0/rho*Qy[4];
  BgradQ[6] = -1.0/rho*Qy[3];
 
  BgradQ[7] = -lam*Qy[1];
  BgradQ[8] = -(lam+2.0*mu)*Qy[1];
 
  BgradQ[9] = -mu*Qy[0];
}

    
void Elastic::MyElasticSolver::coefficientMatrix(const double* const Q,const int normalNonZero,double* Bn) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  // @todo Please implement/augment if required
  Bn[0] = 0.0;
  Bn[1] = 0.0;
  Bn[2] = 0.0;
  Bn[3] = 0.0;
  Bn[4] = 0.0;
  Bn[5] = 0.0;
  Bn[6] = 0.0;
  Bn[7] = 0.0;
  Bn[8] = 0.0;
  Bn[9] = 0.0;
  Bn[10] = 0.0;
  Bn[11] = 0.0;
  Bn[12] = 0.0;
  Bn[13] = 0.0;
  Bn[14] = 0.0;
  Bn[15] = 0.0;
  Bn[16] = 0.0;
  Bn[17] = 0.0;
  Bn[18] = 0.0;
  Bn[19] = 0.0;
  Bn[20] = 0.0;
  Bn[21] = 0.0;
  Bn[22] = 0.0;
  Bn[23] = 0.0;
  Bn[24] = 0.0;

  double lam;  //= 2.0;
  double mu; // = 1.0;
  double rho; // = 1.0;

 
  double cs = 3.464;  // km/s
  double cp = 6.0;    // km/s

  rho = 2.67;   // gm/cm^3
  
  lam = rho * (cp * cp - 2 * cs * cs);  // lambda (MPa)
  mu= rho * cs * cs; 



  double nv[2] = {0.0};
  nv[normalNonZero] = 1.0;
  
  double B1[9][9];
  double B2[9][9];

  for(int i=0; i<5; i++) {
   for(int j=0; j<5; j++) {
     B1[i][j] = 0.0;
     B2[i][j] = 0.0;
   }
 }

  B1[2][0] = -1.0/rho; 
  B1[4][1] = -1.0/rho;


  B1[0][2] = -(2.0*mu + lam); 
  B1[0][3] = -lam;
  
  B1[1][4] = -mu; 
 
  

  B2[4][0] = -1.0/rho; 
  B2[3][1] = -1.0/rho;
  
  B2[1][2] = -lam;
  B2[1][3] = -(2.0*mu + lam);
  
  B2[0][4] = -mu; 
 



  for(int i=0; i<5; i++) {
   for(int j=0; j<5; j++) {
     
     
      Bn[i*5 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
   }
 }
}

//TODO KD
// tell the user what it is



  return true;
}

//TODO KD
void Elastic::MyElasticSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0) {
  //TODO KD // @todo Please implement/augment if required and set bool function

   double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  //double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));

  f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));


  //f = M0*t/(t0*t0)*std::exp(-t/t0);

  x0[0] = 5;
  x0[1] = 5;

  forceVector[0] = 0.0;
  forceVector[1] = 0.0;
  forceVector[2] = 1.*f;
  forceVector[3] = 1.*f;
  forceVector[4] = 0.0;
}