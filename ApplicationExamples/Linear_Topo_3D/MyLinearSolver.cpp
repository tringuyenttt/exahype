#include "MyLinearSolver.h"

#include "MyLinearSolver_Variables.h"

#include "../../ExaHyPE/kernels/KernelUtils.h"

tarch::logging::Log Linear::MyLinearSolver::_log( "Linear::MyLinearSolver" );


void Linear::MyLinearSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Linear::MyLinearSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void Linear::MyLinearSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  double x_0[3]={0.5 , 0.5 ,0.5};
  double exponent=0;

  for(int i=0 ; i < 3 ; i++){
    exponent = exponent + (x[i]-x_0[i])*(x[i]-x_0[i]);
  }
    
  Q[ 0] = std::exp(-exponent/0.01);
  Q[ 1] = 0.0;
  Q[ 2] = 0.0;
  Q[ 3] = 0.0;   // Material parameters:
  Q[ 4] = 1.0;   //rho
  Q[ 5] = 1.484; //c
  Q[ 6] = 0.0;
  Q[ 7] = 0.0;
  Q[ 8] = 0.0;
  Q[ 9] = 0.0;
  Q[10] = 0.0;
  Q[11] = 0.0;
  Q[12] = 0.0;
  Q[13] = 0.0;
  Q[14] = 0.0;
  Q[15] = 0.0;
  Q[16] = 0.0;
  Q[17] = 0.0;
  Q[18] = 0.0;
  
}

void Linear::MyLinearSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required

  double c = Q[5];
  
  lambda[ 0] = c;
  lambda[ 1] = -c;
  lambda[ 2] = 0;
  lambda[ 3] = 0;
}


void Linear::MyLinearSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required

  double u = Q[1];
  double v = Q[2];
  double w = Q[3];

  
  F[0][ 0] = -u;
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;

  F[1][ 0] = -v;
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;

  F[2][ 0] = -w;
  F[2][ 1] = 0.0;
  F[2][ 2] = 0.0;
  F[2][ 3] = 0.0;
}


void Linear::MyLinearSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {

  double p_x = gradQ[0];  
  double u_x = gradQ[1];
  double v_x = gradQ[2];
  double w_x = gradQ[3];

  double p_y = gradQ[4];  
  double u_y = gradQ[5];
  double v_y = gradQ[6];
  double w_y = gradQ[7];

  double p_z = gradQ[8];  
  double u_z = gradQ[9];
  double v_z = gradQ[10];
  double w_z = gradQ[11];    
  
  //  BgradQ[0]= -u_y;
  BgradQ[0]= 0;  
  BgradQ[1]= -p_x;
  BgradQ[2]= 0;
  BgradQ[3]= 0;

  //  BgradQ[4]=-v_y;
  BgradQ[4]=0;  
  BgradQ[5]=0;
  BgradQ[6]=-p_y;
  BgradQ[7]=0;

  //  BgradQ[8]=-w_z;
  BgradQ[8]=0;  
  BgradQ[9]=0;
  BgradQ[10]=0;
  BgradQ[11]=-p_z; 

}



void Linear::MyLinearSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters

  // @todo Please implement/augment if required
  stateOut[ 0] = 0.0;
  stateOut[ 1] = 0.0;
  stateOut[ 2] = 0.0;
  stateOut[ 3] = 0.0;

  fluxOut[ 0] = 0.0;
  fluxOut[ 1] = 0.0;
  fluxOut[ 2] = 0.0;
  fluxOut[ 3] = 0.0;
}


void Linear::MyLinearSolver::coefficientMatrix(const double* const Q,const int d,double* Bn){
  double n[3]={0.0,0.0,0.0};
  n[d]=1.0;

  double B1[4][4];
  double B2[4][4];
  double B3[4][4];

  kernels::idx2 idx_Bn(4,4);

  for (int i =0; i< 4 ; i++){
    for (int j =0; j< 4 ; j++){
      B1[i][j]=0;
      B2[i][j]=0;
      B3[i][j]=0;	
    }
  }

  B1[0][1] = -1.0;
  B1[1][0] = -1.0;

  B2[0][2] = -1.0;
  B2[2][0] = -1.0;

  B3[0][3] = -1.0;
  B3[3][0] = -1.0;  

  for (int i =0; i< 4 ; i++){
    for (int j =0; j< 4 ; j++){
      Bn[idx_Bn(i,j)] = n[0]*B1[i][j]+n[1]*B2[i][j]+n[2]*B3[i][j];
    }
  }
}



void Linear::MyLinearSolver::algebraicSource(const double* const Q,double* S) {

  S[0] = 10*Q[0];
  S[1] = 10*Q[1];
  S[2] = 10*Q[2];
  S[3] = 10*Q[3];
}


void Linear::MyLinearSolver::multiplyMaterialParameterMatrix(const double* const Q, double* rhs){

  double rho = Q[4];  
  double c   = Q[5];
  double mu  = rho*c*c;

  rhs[0]=   mu * rhs[0];
  rhs[1]=1/rho * rhs[1];
  rhs[2]=1/rho * rhs[2];
  rhs[3]=1/rho * rhs[3];
  
}


exahype::solvers::Solver::RefinementControl Linear::MyLinearSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}
