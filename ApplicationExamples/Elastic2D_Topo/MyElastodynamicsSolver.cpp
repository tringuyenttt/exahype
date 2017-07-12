#include "MyElastodynamicsSolver.h"

#include "MyElastodynamicsSolver_Variables.h"

#include "kernels/KernelUtils.h"

#include "kernels/aderdg/generic/Kernels.h"

#include <algorithm>

#include "CurvilinearTransformation.h"



tarch::logging::Log Elastodynamics::MyElastodynamicsSolver::_log( "Elastodynamics::MyElastodynamicsSolver" );


void Elastodynamics::MyElastodynamicsSolver::init(std::vector<std::string>& cmdlineargs) {
   static tarch::logging::Log _log("MyElastodynamicsSolver::init");
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Elastodynamics::MyElastodynamicsSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  //  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PointWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
    return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PatchWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}

void Elastodynamics::MyElastodynamicsSolver::adjustPatchSolution(
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt,
      double* luh) {

  constexpr int basisSize = MyElastodynamicsSolver::Order+1;
  int num_nodes = basisSize;
  int numberOfData=MyElastodynamicsSolver::NumberOfParameters+MyElastodynamicsSolver::NumberOfVariables;

  kernels::idx3 id_xyz(basisSize,basisSize,numberOfData);
  kernels::idx2 id_xy(basisSize,basisSize);


  int nx = 1/dx[0]*num_nodes;
  int ny = 1/dx[1]*num_nodes;


  double* left_bnd_x = new double[ny];

  double* left_bnd_y = new double[ny];

  double* right_bnd_x = new double[ny];
  double* right_bnd_y = new double[ny];

  double* bottom_bnd_x = new double[nx];
  double* bottom_bnd_y = new double[nx];

  double* top_bnd_x = new double[nx];
  double* top_bnd_y = new double[nx];

  double offset_x=cellCentre[0]-0.5*dx[0];
  double offset_y=cellCentre[1]-0.5*dx[1];

  double width_x=dx[0];
  double width_y=dx[1];   

  
  //std::exit(-1);

  getBoundaryCurves(num_nodes, offset_x,  offset_y, width_x,  width_y ,left_bnd_x,left_bnd_y,right_bnd_x,right_bnd_y,bottom_bnd_x,bottom_bnd_y,top_bnd_x,top_bnd_y);


  double* curvilinear_x = new double[num_nodes*num_nodes];
  double* curvilinear_y = new double[num_nodes*num_nodes];

  int i_m =  (offset_x-0.0)/width_x *num_nodes;
  int j_m =  (offset_y-0.0)/width_y *num_nodes;

  // std::cout<< i_m << std::endl;
  // std::cout<< offset_x << std::endl;
  // std::cout<< width_x << std::endl;

  // //std::exit(-1);

  //i0_m = (i-1)*(ndp-1)+1; i0_p = i*(ndp-1) + 1;
  //j0_m = (j-1)*(ndp-1)+1; j0_p = j*(ndp-1) + 1;
	
  int i_p = i_m + num_nodes;
  int j_p = j_m + num_nodes; 
  
  transFiniteInterpolation(nx, ny, j_m, j_p, i_m, i_p, num_nodes,  left_bnd_x,  left_bnd_y,  right_bnd_x,  right_bnd_y,  bottom_bnd_x,  bottom_bnd_y,  top_bnd_x,  top_bnd_y,  curvilinear_x ,  curvilinear_y );

  
  double* gl_vals_x = new double[num_nodes*num_nodes];
  double* gl_vals_y = new double[num_nodes*num_nodes];
  double* jacobian = new double[num_nodes*num_nodes];


  double* q_x = new double[num_nodes*num_nodes];
  double* q_y = new double[num_nodes*num_nodes];
  
  double* r_x = new double[num_nodes*num_nodes];
  double* r_y = new double[num_nodes*num_nodes];  
  
  metricDerivativesAndJacobian(num_nodes,curvilinear_x,curvilinear_y,gl_vals_x,gl_vals_y,q_x,q_y,r_x,r_y,jacobian, width_x, width_y); 

  for (int i=0; i< num_nodes; i++){
      for (int j=0; j< num_nodes; j++){
	double x= gl_vals_x[id_xy(i,j)];
	double y= gl_vals_y[id_xy(i,j)];

	
	//particle velocities
	luh[id_xyz(i,j,0)] = std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01);
	luh[id_xyz(i,j,1)] = std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.01);
	
	// luh[id_xyz(i,j,0)] = x;
	// luh[id_xyz(i,j,1)] = y;
	
	//luh[id_xyz(j,i,0)] = x*0.1;
	
	//Stresses
	luh[id_xyz(i,j,2)] = 0;
	luh[id_xyz(i,j,3)] = 0;
	luh[id_xyz(i,j,4)] = 0;

	luh[id_xyz(i,j,5)] = 2.6;     // gm/cm^3
	luh[id_xyz(i,j,6)] = 2.0;      // km/s
	luh[id_xyz(i,j,7)] = 4.0;      // km/s
	
	if (x > 1.0) {
	  luh[id_xyz(i,j,5)] = 2.7;     // gm/cm^3
	  luh[id_xyz(i,j,6)] = 3.464;    // km/s
	  luh[id_xyz(i,j,7)] = 6.0;      // km/s
	}
	
	luh[id_xyz(i,j,8)] = jacobian[id_xy(i,j)];
	//std::cout << jacobian[id_xy(j,i)] << std::endl;

	luh[id_xyz(i,j,9)] = q_x[id_xy(i,j)];
	luh[id_xyz(i,j,10)] = q_y[id_xy(i,j)];
	luh[id_xyz(i,j,11)] = r_x[id_xy(i,j)];
	luh[id_xyz(i,j,12)] = r_y[id_xy(i,j)];

	// std::cout << q_x[id_xy(j,i)] << std::endl;
	// std::cout << q_y[id_xy(j,i)] << std::endl;
	// std::cout << r_x[id_xy(j,i)] << std::endl;	
	// std::cout << r_y[id_xy(j,i)] << std::endl;
	// std::cout <<  std::endl;
	
	
	luh[id_xyz(i,j,13)] = gl_vals_x[id_xy(i,j)];
	luh[id_xyz(i,j,14)] = gl_vals_y[id_xy(i,j)];	

      }
  }
  
}


void Elastodynamics::MyElastodynamicsSolver::adjustPointSolution(const double* const x,const double w,const double t,const double dt,double* Q) {
  // Dimensions             = 2
  // Number of variables    = 15 + #parameters
  
  // @todo Please implement/augment if required
  // State variables:
  // Q[ 0] = 0.0;
  // Q[ 1] = 0.0;
  // Q[ 2] = 0.0;
  // Q[ 3] = 0.0;
  // Q[ 4] = 0.0;  // Material parameters:
  // Q[ 5] = 0.0;
  // Q[ 6] = 0.0;
  // Q[ 7] = 0.0;
  // Q[ 8] = 0.0;
  // Q[ 9] = 0.0;
  // Q[10] = 0.0;
  // Q[11] = 0.0;
  // Q[12] = 0.0;
  // Q[13] = 0.0;
  // Q[14] = 0.0;
}

void Elastodynamics::MyElastodynamicsSolver::eigenvalues(const double* const Q,const int d,double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 15 + #parameters
  
  // @todo Please implement/augment if required
  double rho  = Q[5];   // km/s
  double cs   = Q[6];   // km/s
  double cp   = Q[7];   // km/s

  //std::cout << rho << cs << cp <<'\n';

  //std::exit(-1);

  //
  lambda[0] = cp;
  lambda[1] = cs;
  lambda[2] = -cp;
  lambda[3] = -cs;
  lambda[4] = 0.0;
 
}


void Elastodynamics::MyElastodynamicsSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 2
  // Number of variables    = 15 + #parameters
  
  // @todo Please implement/augment if required


  double u=Q[0];
  double v=Q[1];
  double sigma_xx=Q[2];
  double sigma_yy=Q[3];
  double sigma_xy=Q[4];

  double jacobian=Q[8];
  double q_x=Q[9];
  double q_y=Q[10];
  double r_x=Q[11];    
  double r_y=Q[12];        
  
  F[0][ 0] = -jacobian*(q_x*sigma_xx+q_y*sigma_xy);
  F[0][ 1] = -jacobian*(q_x*sigma_xy+q_y*sigma_yy);
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;
  F[0][ 4] = 0.0;

  F[1][ 0] = -jacobian*(r_x*sigma_xx+r_y*sigma_xy);
  F[1][ 1] = -jacobian*(r_x*sigma_xy+r_y*sigma_yy);
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;
  F[1][ 4] = 0.0;
}


void Elastodynamics::MyElastodynamicsSolver::boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int normalNonZero,
  const double * const fluxIn,const double* const stateIn,
  double *fluxOut,double* stateOut) {
  // Dimensions             = 2
  // Number of variables    = 15 + #parameters

  // @todo Please implement/augment if required

  double rho  = stateIn[5];   // km/s
  double cs   = stateIn[6];   // km/s
  double cp   = stateIn[7];   // km/s

  


  double n[2] = {0,0};
  n[normalNonZero] = 1.;

  // extract local s-wave and p-wave impedances
  double zp=rho*cp;
  double zs=rho*cs;
  
  stateOut[0] = stateIn[0];
  stateOut[1] = stateIn[1];
  stateOut[2] = stateIn[2];
  stateOut[3] = stateIn[3];
  stateOut[4] = stateIn[4];

  stateOut[5] = stateIn[5];
  stateOut[6] = stateIn[6];
  stateOut[7] = stateIn[7];
  stateOut[8] = stateIn[8];
  stateOut[9] = stateIn[9];

  stateOut[10] = stateIn[10];
  stateOut[11] = stateIn[11];
  stateOut[12] = stateIn[12];
  stateOut[13] = stateIn[13];
  stateOut[14] = stateIn[14];

  fluxOut[0] = fluxIn[0];
  fluxOut[1] = fluxIn[1];
  fluxOut[2] = fluxIn[2];
  fluxOut[3] = fluxIn[3];
  fluxOut[4] = fluxIn[4];

  fluxOut[5] = fluxIn[5];
  fluxOut[6] = fluxIn[6];
  fluxOut[7] = fluxIn[7];
  fluxOut[8] = fluxIn[8];
  fluxOut[9] = fluxIn[9];

  fluxOut[10] = fluxIn[10];
  fluxOut[11] = fluxIn[11];
  fluxOut[12] = fluxIn[12];
  fluxOut[13] = fluxIn[13];
  fluxOut[14] = fluxIn[14];
  
}


exahype::solvers::Solver::RefinementControl Elastodynamics::MyElastodynamicsSolver::refinementCriterion(const double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ){
  kernels::idx2 idx(DIMENSIONS, NumberOfVariables);

  static tarch::logging::Log _log("Elastodynamics::MyElastodynamicsSolver::nonConservativeProduct");
  ReadOnlyVariables vars(Q);

  const double* Qx = &gradQ[0]; // (:,1)
  const double* Qy = &gradQ[5]; // (:,2)

  const double rho  = Q[5];   // km/s
  const double cs   = Q[6];   // km/s
  const double cp   = Q[7];   // km/s

  double jacobian=Q[8];
  double q_x=Q[9];
  double q_y=Q[10];
  double r_x=Q[11];    
  double r_y=Q[12];        
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;

  double u_q=gradQ[0];    
  double v_q=gradQ[1];
  
  double u_r=gradQ[5];
  double v_r=gradQ[6];
  
  BgradQ[0] = 0;
  BgradQ[1] = 0;
  BgradQ[2] = -q_x*u_q;
  BgradQ[3] = -q_y*v_q;
  BgradQ[4] = -(q_y*u_q+q_x*v_q);
  
  BgradQ[5] = 0.0;
  BgradQ[6] = 0.0;
  BgradQ[7] = -r_x*u_r;
  BgradQ[8] = -r_y*v_r;
  BgradQ[9] = -(r_y*u_r+r_x*v_r);


    
  //  BgradQ[0] = -1/rho*Qx[2];
  //  BgradQ[1] = -1/rho*Qx[4];
  // BgradQ[2] = -(lam + 2.0*mu)*Qx[0];
  // BgradQ[3] = -0.0;
  // BgradQ[4] = -mu*Qx[1];
  // BgradQ[5] = -1/rho*Qy[4];
  // BgradQ[6] = -1/rho*Qy[3];
  // BgradQ[7] = -0.0;
  // BgradQ[8] = -(lam + 2*mu)*Qy[1];
  // BgradQ[9] = -mu*Qy[0];



}


void Elastodynamics::MyElastodynamicsSolver::coefficientMatrix(const double* const Q,const int d,double* Bn){
  static tarch::logging::Log _log("MyElastodynamicsSolver::coefficientMatrix");

  const double rho  = Q[5];   // km/s
  const double cs   = Q[6];   // km/s
  const double cp   = Q[7];   // km/s
    
  double mu = rho*cs*cs;
  double lam = rho*cp*cp - 2*mu;
  
  double nv[2] = {0.0, 0.0};

  nv[d] = 1.0;
  
  double B1[5][5];
  double B2[5][5];
   
  for(int i=0; i<5; i++) {
    for(int j=0; j<5; j++) {
      B1[i][j] = 0.0;
      B2[i][j] = 0.0;
    }
  }
  
  B1[0][2] = -1/rho; 
  B1[1][4] = -1/rho;
  B1[2][0] = -(2*mu + lam); 
  B1[4][1] = -mu;

  B2[0][4] = -1/rho; 
  B2[1][3] = -1/rho;
  B2[3][1] = -(2*mu + lam); 
  B2[4][0] = -mu;
  
  for(int i=0; i<5; i++) {
    for(int j=0; j<5; j++) {
      
      Bn[i*5 + j] = nv[0]*B1[i][j] + nv[1]*B2[i][j];
    }
  }

}


void Elastodynamics::MyElastodynamicsSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex,bool isBounaryFace, int faceIndex){

  constexpr int numberOfVariables  = MyElastodynamicsSolver::NumberOfVariables;
  constexpr int numberOfVariables2 = numberOfVariables*numberOfVariables;
  constexpr int numberOfParameters = MyElastodynamicsSolver::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int basisSize          = MyElastodynamicsSolver::Order+1;
  constexpr int order              = basisSize - 1;
  
  kernels::idx2 idx_QLR(basisSize, numberOfData);

  kernels::idx2 idx_FLR(basisSize, NumberOfVariables);

  double n[2]={0,0};
  n[normalNonZeroIndex]=1;
  
  for (int i = 0; i < basisSize; i++) {
    double qm_x=QL[idx_QLR(i,9)];
    double qm_y=QL[idx_QLR(i,10)];
    double rm_x=QL[idx_QLR(i,11)];
    double rm_y=QL[idx_QLR(i,12)];

    double qp_x=QR[idx_QLR(i,9)];
    double qp_y=QR[idx_QLR(i,10)];
    double rp_x=QR[idx_QLR(i,11)];
    double rp_y=QR[idx_QLR(i,12)];

    double n_p[3]={0,0,0};
    double n_m[3]={0,0,0};
    
    double m_p[3]={0,0,0};
    double m_m[3]={0,0,0};

    double l_p[3]={0,0,0};
    double l_m[3]={0,0,0};
    
    
    double norm_p_qr;
    double norm_m_qr;
    
    if (normalNonZeroIndex == 0){
      
      norm_m_qr = std::sqrt(qm_x*qm_x + qm_y*qm_y);
      n_m[0] = qm_x/norm_m_qr;
      n_m[1] = qm_y/norm_m_qr;

      norm_p_qr = std::sqrt(qp_x*qp_x + qp_y*qp_y);
      n_p[0] = qp_x/norm_p_qr;
      n_p[1] = qp_y/norm_p_qr;
      
      m_m[0] = n_m[1];
      m_m[1] =-n_m[0];
      
      m_p[0] = n_p[1];
      m_p[1] =-n_p[0];
      }
      
    if (normalNonZeroIndex == 1){

      norm_m_qr = std::sqrt(rm_x*rm_x + rm_y*rm_y);
      n_m[0] = rm_x/norm_m_qr;
      n_m[1] = rm_y/norm_m_qr;

      norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y);
      n_p[0] = rp_x/norm_p_qr;
      n_p[1] = rp_y/norm_p_qr;

      m_m[0] = n_m[1];
      m_m[1] =-n_m[0];

      m_p[0] = n_p[1];
      m_p[1] =-n_p[0];
	
      // norm_qr = std::sqrt(r_x*r_x + r_y*r_y);
      // n[0] = r_x/norm_qr;
      // n[1] = r_y/norm_qr;
    }

    double rho_m  = QL[idx_QLR(i,5)];   // km/s
    double cs_m   = QL[idx_QLR(i,6)];   // km/s
    double cp_m   = QL[idx_QLR(i,7)];   // km/s

    double rho_p  = QR[idx_QLR(i,5)];   // km/s
    double cs_p   = QR[idx_QLR(i,6)];   // km/s
    double cp_p   = QR[idx_QLR(i,7)];   // km/s
    
    //std::cout << rho_p << cs_p << cp_p <<numberOfData<<'\n';

    //std::exit(-1);
    
    double mu_m = rho_m*cs_m*cs_m;
    double lam_m = rho_m*cp_m*cp_m - 2*mu_m;

    double mu_p = rho_p*cs_p*cs_p;
    double lam_p = rho_p*cp_p*cp_p - 2*mu_p;
  

    // extract tractions and particle velocities
    double sigma_m_xx =  QL[idx_QLR(i,2)];
    double sigma_m_yy =  QL[idx_QLR(i,3)];
    double sigma_m_xy =  QL[idx_QLR(i,4)];

    double sigma_p_xx =  QR[idx_QLR(i,2)];
    double sigma_p_yy =  QR[idx_QLR(i,3)];
    double sigma_p_xy =  QR[idx_QLR(i,4)];

    double Tx_m = n_m[0]*sigma_m_xx + n_m[1]*sigma_m_xy;
    double Ty_m = n_m[0]*sigma_m_xy + n_m[1]*sigma_m_yy;

    double Tx_p = n_p[0]*sigma_p_xx + n_p[1]*sigma_p_xy;
    double Ty_p = n_p[0]*sigma_p_xy + n_p[1]*sigma_p_yy;

    double vx_m = QL[idx_QLR(i,0)];
    double vy_m = QL[idx_QLR(i,1)];

    double vx_p = QR[idx_QLR(i,0)];
    double vy_p = QR[idx_QLR(i,1)];

    localBasis(n_m, m_m, l_m, 2);
    localBasis(n_p, m_p, l_p, 2);    

    // rotate tractions and particle velocities into orthogonal coordinates: n, m
    double Tn_m= Tx_m*n_m[0] + Ty_m*n_m[1];
    double Tm_m= Tx_m*m_m[0] + Ty_m*m_m[1];
    
    double Tn_p= Tx_p*n_p[0] + Ty_p*n_p[1];
    double Tm_p= Tx_p*m_p[0] + Ty_p*m_p[1];
    
    double vn_m= vx_m*n_m[0] + vy_m*n_m[1];
    double vm_m= vx_m*m_m[0] + vy_m*m_m[1];
    
    double vn_p= vx_p*n_p[0] + vy_p*n_p[1];
    double vm_p= vx_p*m_p[0] + vy_p*m_p[1];

   
    // extract local s-wave and p-wave impedances
    double zs_p=rho_p*cs_p;
    double zs_m=rho_m*cs_m;

    double zp_p=rho_p*cp_p;
    double zp_m=rho_m*cp_m;

    // impedance must be greater than zero !
    if (zs_p <= 0.0 || zs_m <= 0.0 || zp_p <= 0.0 || zp_m <= 0.0){
      std::cout<<zs_p<<' '<<zs_m<<' '<<zp_p<<' '<<zp_m<<'\n';
      std::cout<<' Impedance must be greater than zero ! '<<'\n';
      std::exit(-1);
    }

    
    // generate interface data preserving the amplitude of the outgoing charactertritics
    // and satisfying interface conditions exactly.
    double vn_hat_p=0;
    double vm_hat_p=0;

    double vn_hat_m=0;
    double vm_hat_m=0;

    double Tn_hat_p=0;
    double Tm_hat_p=0;

    double Tn_hat_m=0;
    double Tm_hat_m=0;

    // data is generated by solving a local Riemann problem and contraining the solutions against
    // physical interface conditions
    

    if (isBounaryFace) {
      if (faceIndex == 0) {
	
	double r = 0.;
	
	riemannSolver_BC0(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
	riemannSolver_BC0(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
	
	riemannSolver_BC0(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
	riemannSolver_BC0(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
	
      }
      
      
      if (faceIndex == 1) {
	
	
	double r = 0.;
	
	riemannSolver_BCn(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
	riemannSolver_BCn(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
	
	riemannSolver_BCn(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
	riemannSolver_BCn(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
        
      }
      
      
      if (faceIndex == 2) {
	
	double r = 0.;
	
	riemannSolver_BC0(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
	riemannSolver_BC0(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
	
	riemannSolver_BC0(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
	riemannSolver_BC0(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
	
	
      }
      
      if (faceIndex == 3) {
	
	double r = 1.;
	
	riemannSolver_BCn(vn_p, Tn_p, zp_p, r, vn_hat_p, Tn_hat_p);
	riemannSolver_BCn(vm_p, Tm_p, zs_p, r, vm_hat_p, Tm_hat_p);
	
	riemannSolver_BCn(vn_m, Tn_m, zp_m, r, vn_hat_m, Tn_hat_m);
	riemannSolver_BCn(vm_m, Tm_m, zs_m, r, vm_hat_m, Tm_hat_m);
        
      }
    }
    else {
      
      riemannSolver_Nodal(vn_p,vn_m, Tn_p, Tn_m, zp_p , zp_m, vn_hat_p , vn_hat_m, Tn_hat_p, Tn_hat_m);
      riemannSolver_Nodal(vm_p,vm_m, Tm_p, Tm_m, zs_p , zs_m, vm_hat_p , vm_hat_m, Tm_hat_p, Tm_hat_m);
    }
    // generate fluctuations in the local basis coordinates: n, m
    double FLn = 0.5*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FLm = 0.5*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));

    double FRn = 0.5*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FRm = 0.5*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));

    double FL_n = 0.5/zp_m*(zp_m*(vn_m-vn_hat_m) + (Tn_m-Tn_hat_m));
    double FL_m = 0.5/zs_m*(zs_m*(vm_m-vm_hat_m) + (Tm_m-Tm_hat_m));

    double FR_n = 0.5/zp_p*(zp_p*(vn_p-vn_hat_p) - (Tn_p-Tn_hat_p));
    double FR_m = 0.5/zs_p*(zs_p*(vm_p-vm_hat_p) - (Tm_p-Tm_hat_p));

    // rotate back to the physical coordinates x, y
    double FLx = n_m[0]*FLn + m_m[0]*FLm;
    double FLy = n_m[1]*FLn + m_m[1]*FLm;

    double FRx = n_p[0]*FRn + m_p[0]*FRm;
    double FRy = n_p[1]*FRn + m_p[1]*FRm;

    double FL_x = n_m[0]*FL_n + m_m[0]*FL_m;
    double FL_y = n_m[1]*FL_n + m_m[1]*FL_m;

    double FR_x = n_p[0]*FR_n + m_p[0]*FR_m;
    double FR_y = n_p[1]*FR_n + m_p[1]*FR_m;
     
    // construct flux fluctuation vectors obeying the eigen structure of the PDE
    // and choose physically motivated penalties such that we can prove
    // numerical stability.

    FR[idx_FLR(i, 0)] = norm_p_qr/rho_p*FRx;
    FL[idx_FLR(i, 0)] = norm_m_qr/rho_m*FLx;
    
    FR[idx_FLR(i, 1)] = norm_p_qr/rho_p*FRy;
    FL[idx_FLR(i, 1)] = norm_m_qr/rho_m*FLy;

    FL[idx_FLR(i, 2)] = norm_m_qr*((2*mu_m+lam_m)*n_m[0]*FL_x + lam_m*n_m[1]*FL_y);
    FL[idx_FLR(i, 3)] = norm_m_qr*((2*mu_m+lam_m)*n_m[1]*FL_y + lam_m*n_m[0]*FL_x);

    FR[idx_FLR(i, 2)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[0]*FR_x + lam_p*n_p[1]*FR_y);
    FR[idx_FLR(i, 3)] = -norm_p_qr*((2*mu_p+lam_p)*n_p[1]*FR_y + lam_p*n_p[0]*FR_x);


    FL[idx_FLR(i, 4)] =  norm_m_qr*mu_m*(n_m[1]*FL_x + n_m[0]*FL_y);
    FR[idx_FLR(i, 4)] = -norm_p_qr*mu_p*(n_p[1]*FR_x + n_p[0]*FR_y);
    
  }
  
}


void Elastodynamics::MyElastodynamicsSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0){
  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.1;
  double f = 0.0;
  double M0 = 1000.0;
  

  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-pow((t-t0)/(std::sqrt(2)*sigma),2.0)));
  
  //f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
  
  f = M0*t/(t0*t0)*std::exp(-t/t0);
  
  x0[0] = 2.0;
  x0[1] = 15.0;
  
  forceVector[0] = 0.0;
  forceVector[1] = 0.0;
  forceVector[2] = 0;
  forceVector[3] = 0;
  forceVector[4] = f;

}


void Elastodynamics::MyElastodynamicsSolver::riemannSolver_Nodal(double v_p,double v_m, double sigma_p, double sigma_m, double z_p , double z_m, double& v_hat_p , double& v_hat_m, double& sigma_hat_p, double& sigma_hat_m){
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


void Elastodynamics::MyElastodynamicsSolver::Gram_Schmidt(double* y, double* z){
  //Gram Schmidt orthonormalization
 
  double  a_yz = y[0]*z[0] + y[1]*z[1] + y[2]*z[2];

  for (int i = 0; i< 3; i++){
    z[i] = z[i] - a_yz*y[i];
    z[i] = 1.0/std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2])*z[i];
  }

}

void Elastodynamics::MyElastodynamicsSolver::localBasis(double* n, double * m, double* l, int d){

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


void Elastodynamics::MyElastodynamicsSolver::multiplyMaterialParameterMatrix(const double* const Q, double* rhs){

  // double jacobian= Q[4];

  const double rho  = Q[5];   // km/s
  const double cs   = Q[6];   // km/s
  const double cp   = Q[7];   // km/s
  double jacobian = Q[8];
  
  double mu = rho*cs*cs;
  double lam = rho*cp*cp-2.0*mu;


  rhs[0] = 1/(rho*jacobian)*rhs[0];
  rhs[1] = 1/(rho*jacobian)*rhs[1];
  
  double rhs_2= (2*mu+lam)*rhs[2]+lam*rhs[3];
  double rhs_3= (2*mu+lam)*rhs[3]+lam*rhs[2];
  
  rhs[2]=rhs_2;
  rhs[3]=rhs_3;  
  rhs[4]=mu*rhs[4];


  rhs[5] = 1/(rho*jacobian)*rhs[5];
  rhs[6] = 1/(rho*jacobian)*rhs[6];
  
  double rhs_7= (2*mu+lam)*rhs[7]+lam*rhs[8];
  double rhs_8= (2*mu+lam)*rhs[8]+lam*rhs[7];
  
  rhs[7]=rhs_7;
  rhs[8]=rhs_8;  
  rhs[9]=mu*rhs[9];

}

void Elastodynamics::MyElastodynamicsSolver::riemannSolver_BC0(double v, double sigma, double z,  double r, double& v_hat, double& sigma_hat){
  
   double p = 0.5*(z*v + sigma);

   v_hat = (1+r)/z*p;
   sigma_hat = (1-r)*p;

 }


void Elastodynamics::MyElastodynamicsSolver::riemannSolver_BCn(double v,double sigma, double z, double r, double& v_hat, double& sigma_hat){
  
   double q = 0.5*(z*v - sigma);

   v_hat = (1+r)/z*q;
   sigma_hat = -(1-r)*q;

 }
