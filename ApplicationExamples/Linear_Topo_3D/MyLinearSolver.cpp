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
  
  //BgradQ[0]= -u_x;
  BgradQ[0]= 0;  
  BgradQ[1]= -p_x;
  BgradQ[2]= 0;
  BgradQ[3]= 0;

  //BgradQ[4]=-v_y;
  BgradQ[4]=0;  
  BgradQ[5]=0;
  BgradQ[6]=-p_y;
  BgradQ[7]=0;

  //BgradQ[8]=-w_z;
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
  
  for(int i = 0 ; i< 19; i++){
    stateOut[ i] = stateIn[i];
  }


  fluxOut[ 0] = fluxIn[ 0];
  fluxOut[ 1] = fluxIn[ 1];
  fluxOut[ 2] = fluxIn[ 2];
  fluxOut[ 3] = fluxIn[ 3];  

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

  rhs[4]=   mu * rhs[4];
  rhs[5]=1/rho * rhs[5];
  rhs[6]=1/rho * rhs[6];
  rhs[7]=1/rho * rhs[7];

  rhs[8]=   mu * rhs[8];
  rhs[9]=1/rho * rhs[9];
  rhs[10]=1/rho * rhs[10];
  rhs[11]=1/rho * rhs[11];
  
}


void Linear::MyLinearSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex,bool isBoundaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyLinearSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyLinearSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyLinearSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx3 idx_QLR(basisSize,basisSize,numberOfData);

  kernels::idx3 idx_FLR(basisSize,basisSize,numberOfVariables);

  double n[3]={0,0,0};

  n[normalNonZeroIndex]=1;


  double n_p[3]={0,0,0};
  double n_m[3]={0,0,0};

  double m_p[3]={0,0,0};
  double m_m[3]={0,0,0};

  double l_p[3]={0,0,0};  
  double l_m[3]={0,0,0};

  double norm_p_qr=1.0;
  double norm_m_qr=1.0;

  //std::cout<<isBoundaryFace<<std::endl;
  

  //std::cout<<isBoundaryFace<<std::endl;

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

      double c_p=QR[idx_QLR(i,j,5)];
      double rho_p=QR[idx_QLR(i,j,4)];
      double mu_p=c_p*c_p*rho_p;

      double c_m=QL[idx_QLR(i,j,5)];
      double rho_m=QL[idx_QLR(i,j,4)];
      double mu_m=c_m*c_m*rho_m;

      // double qm_x=QL[idx_QLR(i,4)];
      // double qm_y=QL[idx_QLR(i,5)];
      // double rm_x=QL[idx_QLR(i,6)];
      // double rm_y=QL[idx_QLR(i,7)];

      // double qp_x=QR[idx_QLR(i,4)];
      // double qp_y=QR[idx_QLR(i,5)];
      // double rp_x=QR[idx_QLR(i,6)];
      // double rp_y=QR[idx_QLR(i,7)];
      
      // if (normalNonZeroIndex == 0){
	
      // 	norm_m_qr = std::sqrt(qm_x*qm_x + qm_y*qm_y);
      // 	n_m[0] = qm_x/norm_m_qr;
      // 	n_m[1] = qm_y/norm_m_qr;

      // 	norm_p_qr = std::sqrt(qp_x*qp_x + qp_y*qp_y);
      // 	n_p[0] = qp_x/norm_p_qr;
      // 	n_p[1] = qp_y/norm_p_qr;

      // 	 m_m[0] = n_m[1];
      // 	 m_m[1] =-n_m[0];

      // 	 m_p[0] = n_p[1];
      // 	 m_p[1] =-n_p[0];

      // }
      
      // if (normalNonZeroIndex == 1){

      // 	norm_m_qr = std::sqrt(rm_x*rm_x + rm_y*rm_y);
      // 	n_m[0] = rm_x/norm_m_qr;
      // 	n_m[1] = rm_y/norm_m_qr;

      // 	norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y);
      // 	n_p[0] = rp_x/norm_p_qr;
      // 	n_p[1] = rp_y/norm_p_qr;

      // 	m_m[0] = n_m[1];
      // 	m_m[1] =-n_m[0];

      // 	m_p[0] = n_p[1];
      // 	m_p[1] =-n_p[0];
	
      // 	// norm_qr = std::sqrt(r_x*r_x + r_y*r_y);
      // 	// n[0] = r_x/norm_qr;
      // 	// n[1] = r_y/norm_qr;
      // }

      for(int k = 0 ; k< 3 ;k++){
	n_m[k] = n[k];
	n_p[k] = n[k];
      }

      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);     
	

      // 	norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y);
      // 	n_p[0] = rp_x/norm_p_qr;
      // 	n_p[1] = rp_y/norm_p_qr;

      // 	m_m[0] = n_m[1];
      // 	m_m[1] =-n_m[0];

      // 	m_p[0] = n_p[1];
      // 	m_p[1] =-n_p[0];

      
      
      double v_m=QL[idx_QLR(i,j,1)]*n_m[0]+QL[idx_QLR(i,j,2)]*n_m[1]+QL[idx_QLR(i,j,3)]*n_m[2];
      double v_p=QR[idx_QLR(i,j,1)]*n_p[0]+QR[idx_QLR(i,j,2)]*n_p[1]+QR[idx_QLR(i,j,3)]*n_p[2];

      double sigma_m = QL[idx_QLR(i,j,0)];
      double sigma_p = QR[idx_QLR(i,j,0)];

      double z_p=rho_p*c_p;
      double z_m=rho_m*c_m;

      double v_hat_p=0;
      double v_hat_m=0;
      double sigma_hat_p=0;
      double sigma_hat_m=0;

      if (isBoundaryFace) {
	if (faceIndex == 0) {
	
	  double r = 0.;
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	}
      
      
	if (faceIndex == 1) {
	  double r = 0.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	}
      
      
	if (faceIndex == 2) {
	
	  double r = 0.;
	
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	
	}
      
	if (faceIndex == 3) {
	
	  double r = 1.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
        
	}

	if (faceIndex == 4) {
	
	  double r = 0.;
	
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	
	}
      
	if (faceIndex == 5) {
	
	  double r = 0.;
	
	  riemannSolver_BCn(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BCn(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
        
	}


	
      }
      else {
	riemannSolver_Nodal(v_p,v_m, sigma_p, sigma_m, z_p , z_m, v_hat_p , v_hat_m, sigma_hat_p, sigma_hat_m);      
      }
      




      FR[idx_FLR(i,j,0)] = -norm_p_qr*0.5*mu_p*((v_p-v_hat_p) - (sigma_p-sigma_hat_p)/z_p);
      FL[idx_FLR(i,j,0)] =  norm_m_qr*0.5*mu_m*((v_m-v_hat_m) + (sigma_m-sigma_hat_m)/z_m);

      FR[idx_FLR(i,j,1)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[0];
      FL[idx_FLR(i,j,1)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[0];

      FR[idx_FLR(i,j,2)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[1];
      FL[idx_FLR(i,j,2)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[1];

      FR[idx_FLR(i,j,3)] = norm_p_qr*0.5/rho_p*(z_p*(v_p-v_hat_p) - (sigma_p-sigma_hat_p))*n_p[2];
      FL[idx_FLR(i,j,3)] = norm_m_qr*0.5/rho_m*(z_m*(v_m-v_hat_m) + (sigma_m-sigma_hat_m))*n_m[2];
      

    }
  }
}



void Linear::MyLinearSolver::Gram_Schmidt(double* y, double* z){
  //Gram Schmidt orthonormalization
 
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
    z[i] = 1.0/std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2])*z[i];
  }

}

void Linear::MyLinearSolver::localBasis(double* n, double * m, double* l, int d){

  if (d == 2)
    {
      l[0] = 0.;
      l[1] = 0.;
      l[2] = 1.0;
      
      m[0] = n[1]*l[2]-n[2]*l[1];
      m[1] = -(n[0]*l[2]-n[2]*l[0]);
      m[2] = n[0]*l[1]-n[1]*l[0];
      
    }
  
  
  if (d == 3)
    {
      double tol, diff_norm1, diff_norm2;
      
      tol = 1e-12;
      m[0] = 0.;
      m[1] = 1.;
      m[2] = 0.;
      
      diff_norm1 =  std::sqrt(pow(n[0]-m[0],2) + pow(n[1]-m[1], 2) + pow(n[2]-m[2], 2));
      
      diff_norm2 =  std::sqrt(pow(n[0]+m[0],2) + pow(n[1]+m[1], 2) + pow(n[2]+m[2], 2));
      
      if (diff_norm1 >= tol && diff_norm2 >= tol){
	Gram_Schmidt(n, m);}	else
	{
	  m[0] = 0.;
	  m[1] = 0.;
	  m[2] = 1.;
	  
	  Gram_Schmidt(n, m);
	  
	}
      
      l[0] = n[1]*m[2]-n[2]*m[1];
      l[1] = -(n[0]*m[2]-n[2]*m[0]);
      l[2] = n[0]*m[1]-n[1]*m[0];
      
    }
  
}



void Linear::MyLinearSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
   double p=0;
   double q=0;
   double phi=0;
   double v_hat=0;
   double eta=0;

   p=z_m*v_p + sigma_p;
   q=z_p*v_m - sigma_m;

   eta=(z_p*z_m)/(z_p+z_m);

   phi= eta*(p/z_p - q/z_m);

   sigma_hat_p=phi;
   sigma_hat_m=phi;

   v_hat_m=(p-phi)/z_p;
   v_hat_p=(q+phi)/z_m;

 }


void Linear::MyLinearSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

}


void Linear::MyLinearSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

}



exahype::solvers::Solver::RefinementControl Linear::MyLinearSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}
