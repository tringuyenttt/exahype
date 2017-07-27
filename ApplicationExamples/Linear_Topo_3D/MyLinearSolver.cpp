#include "MyLinearSolver.h"

#include "MyLinearSolver_Variables.h"

#include "../../ExaHyPE/kernels/KernelUtils.h"
#include "../../ExaHyPE/kernels/DGMatrices.h"

#include "CurvilinearTransformation.h"

tarch::logging::Log Linear::MyLinearSolver::_log( "Linear::MyLinearSolver" );


void Linear::MyLinearSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

exahype::solvers::ADERDGSolver::AdjustSolutionValue Linear::MyLinearSolver::useAdjustSolution(const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,const double t,const double dt) const {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0) ? exahype::solvers::ADERDGSolver::AdjustSolutionValue::PatchWisely : exahype::solvers::ADERDGSolver::AdjustSolutionValue::No;
}



void Linear::MyLinearSolver::adjustPatchSolution(
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt,
      double* luh) {

  constexpr int basisSize = MyLinearSolver::Order+1;
  int num_nodes = basisSize;
  int numberOfData=MyLinearSolver::NumberOfParameters+MyLinearSolver::NumberOfVariables;

  int nx = std::ceil((1/dx[0]))*(num_nodes-1)+1;
  int ny = std::ceil((1/dx[1]))*(num_nodes-1)+1;
  int nz = std::ceil((1/dx[2]))*(num_nodes-1)+1;

  
 // std::cout << nx << " "<< ny << " "<< nz <<std::endl;

 // std::exit(-1);

  kernels::idx4 id_4(basisSize,basisSize,basisSize,numberOfData);
  kernels::idx3 id_3(basisSize,basisSize,basisSize);

  double* left_bnd_x = new double[ny*nz];
  double* left_bnd_y = new double[ny*nz];
  double* left_bnd_z = new double[ny*nz];  

  double* right_bnd_x = new double[ny*nz];
  double* right_bnd_y = new double[ny*nz];
  double* right_bnd_z = new double[ny*nz];  

  double* bottom_bnd_x = new double[nx*nz];
  double* bottom_bnd_y = new double[nx*nz];
  double* bottom_bnd_z = new double[nx*nz];  

  double* top_bnd_x = new double[nx*nz];
  double* top_bnd_y = new double[nx*nz];
  double* top_bnd_z = new double[nx*nz];  

  double* front_bnd_x = new double[nx*ny];
  double* front_bnd_y = new double[nx*ny];
  double* front_bnd_z = new double[nx*ny];  

  double* back_bnd_x = new double[nx*ny];
  double* back_bnd_y = new double[nx*ny];
  double* back_bnd_z = new double[nx*ny];  
  

  double offset_x=cellCentre[0]-0.5*dx[0];
  double offset_y=cellCentre[1]-0.5*dx[1];
  double offset_z=cellCentre[2]-0.5*dx[2];  

  double width_x=dx[0];
  double width_y=dx[1];
  double width_z=dx[2];
  
  // getBoundaryCurves3D( num_nodes,
  // 		       offset_x,  offset_y,  offset_z,
  // 		       width_x,  width_y ,  width_z ,
  // 		       left_bnd_x,  left_bnd_y,  left_bnd_z,
  // 		       right_bnd_x,  right_bnd_y,  right_bnd_z,
  // 		       bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  // 		       top_bnd_x,  top_bnd_y,  top_bnd_z,
  // 		       front_bnd_x,  front_bnd_y,  front_bnd_z,
  // 		       back_bnd_x,  back_bnd_y,  back_bnd_z);

  getBoundaryCurves3D_fixedTopFace( num_nodes,
  				    offset_x,  offset_y,  offset_z,
  		       width_x,  width_y ,  width_z ,
  		       left_bnd_x,  left_bnd_y,  left_bnd_z,
  		       right_bnd_x,  right_bnd_y,  right_bnd_z,
  		       bottom_bnd_x,  bottom_bnd_y,  bottom_bnd_z,
  		       top_bnd_x,  top_bnd_y,  top_bnd_z,
  		       front_bnd_x,  front_bnd_y,  front_bnd_z,
  		       back_bnd_x,  back_bnd_y,  back_bnd_z);


  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); //left right

  double* curvilinear_x = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_y = new double[num_nodes*num_nodes*num_nodes];
  double* curvilinear_z = new double[num_nodes*num_nodes*num_nodes];  

  int i_m;
  int j_m;
  int k_m;


  i_m =  std::round((offset_x/width_x) *(num_nodes-1));
  j_m =  std::round((offset_y/width_y) *(num_nodes-1));
  k_m =  std::round((offset_z/width_z) *(num_nodes-1));

  // std::cout<< width_x<< "  " <<  width_y<< "  " <<  width_z<< "  " << std::endl;
  // std::cout<< offset_x<< "  " << offset_y<< "  " << offset_z<< "  " << std::endl;

  //  std::cout<< std::endl;

  // if (int(offset_x/width_x) == 0)
  //   {
  //     i_m =  (offset_x/width_x);
  //   }

  //  if (int(offset_y/width_y) == 0)
  //   {
  //     j_m =  (offset_y/width_y);
  //   }

  //   if (int(offset_z/width_z) == 0)
  //   {
  //     k_m =  (offset_z/width_z);
  //   }
	
  int i_p = i_m + num_nodes;
  int j_p = j_m + num_nodes;
  int k_p = k_m + num_nodes;   

  
  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,
  			      left_bnd_x,
  			      right_bnd_x,
  			      bottom_bnd_x,
  			      top_bnd_x,
  			      front_bnd_x,
  			      back_bnd_x,
  			      curvilinear_x
  			      );

  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,
  			      left_bnd_y,
  			      right_bnd_y,
  			      bottom_bnd_y,
  			      top_bnd_y,
  			      front_bnd_y,
  			      back_bnd_y,
  			      curvilinear_y
  			      );

  //  double right_bnd_z_block = right_bnd_z[k_m:k_p]
  
  transFiniteInterpolation3D( nx,  ny,  nz,
  			      k_m,  k_p ,
  			      j_m,  j_p ,
  			      i_m,  i_p ,
  			      num_nodes,
			      width_x,width_y,width_z,			      
  			      left_bnd_z,
  			      right_bnd_z,
  			      bottom_bnd_z,
  			      top_bnd_z,
  			      front_bnd_z,
  			      back_bnd_z,
  			      curvilinear_z
  			      );

  
  double* gl_vals_x = new double[num_nodes*num_nodes*num_nodes];
  double* gl_vals_y = new double[num_nodes*num_nodes*num_nodes];
  double* gl_vals_z = new double[num_nodes*num_nodes*num_nodes];

  double* jacobian = new double[num_nodes*num_nodes*num_nodes];

  double* q_x = new double[num_nodes*num_nodes*num_nodes];
  double* q_y = new double[num_nodes*num_nodes*num_nodes];
  double* q_z = new double[num_nodes*num_nodes*num_nodes];
  
  double* r_x = new double[num_nodes*num_nodes*num_nodes];
  double* r_y = new double[num_nodes*num_nodes*num_nodes];
  double* r_z = new double[num_nodes*num_nodes*num_nodes];

  double* s_x = new double[num_nodes*num_nodes*num_nodes];
  double* s_y = new double[num_nodes*num_nodes*num_nodes];
  double* s_z = new double[num_nodes*num_nodes*num_nodes];  

  
  metricDerivativesAndJacobian3D(num_nodes,
  				 curvilinear_x,  curvilinear_y,  curvilinear_z,
  				 gl_vals_x,  gl_vals_y,  gl_vals_z,
  				 q_x,  q_y,  q_z,
  				 r_x,  r_y,  r_z,
  				 s_x,  s_y,  s_z,				  
  				 jacobian,
  				 width_x,  width_y,  width_z
  				 );


  for (int k=0; k< num_nodes; k++){
    for (int j=0; j< num_nodes; j++){
      for (int i=0; i< num_nodes; i++){


	
	double x= gl_vals_x[id_3(k,j,i)];
	double y= gl_vals_y[id_3(k,j,i)];
	double z= gl_vals_z[id_3(k,j,i)];


	std::cout << std::endl;	
	std::cout <<"x" << x - (offset_x+width_x*kernels::gaussLegendreNodes[num_nodes-1][i]) << std::endl;
	std::cout <<"y" << y - (offset_y+width_y*kernels::gaussLegendreNodes[num_nodes-1][j]) << std::endl;
	std::cout <<"z" << z - (offset_z+width_z*kernels::gaussLegendreNodes[num_nodes-1][k]) << std::endl;

	std::cout << std::endl;	

	//Pressure
	//luh[id_4(k,j,i,0)]  = std::exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01);
	luh[id_4(k,j,i,0)] = 0.0;
	
	// //Velocity
	luh[id_4(k,j,i,1)]  = 0;
	luh[id_4(k,j,i,2)]  = 0;
	luh[id_4(k,j,i,3)]  = 0;

	luh[id_4(k,j,i,4)]  = 1.0;   //rho
	luh[id_4(k,j,i,5)]  = 1.484; //c	
	  
	luh[id_4(k,j,i,6)]  = jacobian[id_3(k,j,i)];

	luh[id_4(k,j,i,7)]  = q_x[id_3(k,j,i)];
	luh[id_4(k,j,i,8)]  = q_y[id_3(k,j,i)];
	luh[id_4(k,j,i,9)]  = q_z[id_3(k,j,i)];
	  
	luh[id_4(k,j,i,10)] = r_x[id_3(k,j,i)];
	luh[id_4(k,j,i,11)] = r_y[id_3(k,j,i)];
	luh[id_4(k,j,i,12)] = r_z[id_3(k,j,i)];
	  
	luh[id_4(k,j,i,13)] = s_x[id_3(k,j,i)];
	luh[id_4(k,j,i,14)] = s_y[id_3(k,j,i)];
	luh[id_4(k,j,i,15)] = s_z[id_3(k,j,i)];
	  
	luh[id_4(k,j,i,16)] = gl_vals_x[id_3(k,j,i)];
	luh[id_4(k,j,i,17)] = gl_vals_y[id_3(k,j,i)];
	luh[id_4(k,j,i,18)] = gl_vals_z[id_3(k,j,i)];

      }
    }
  }
  
  // //std::cout << jacobian[id_3(j,i)] << std::endl;
	// // std::cout << q_x[id_xy(j,i)] << std::endl;
	// // std::cout << q_y[id_xy(j,i)] << std::endl;
	// // std::cout << r_x[id_xy(j,i)] << std::endl;	
	// // std::cout << r_y[id_xy(j,i)] << std::endl;
	// // std::cout <<  std::endl;
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
    
  Q[ 0] = 0.0; //std::exp(-exponent/0.01);
  Q[ 1] = 0.0;
  Q[ 2] = 0.0;
  Q[ 3] = 0.0;   // Material parameters:
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

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];
  
  lambda[ 0] = std::sqrt(q_x*q_x + q_y*q_y + q_z*q_z)*c;
  lambda[ 1] = std::sqrt(r_x*r_x + r_y*r_y + r_z*r_z)*c;
  lambda[ 2] = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z)*c;
  lambda[ 3] = 0;
}


void Linear::MyLinearSolver::flux(const double* const Q,double** F) {
  // Dimensions             = 3
  // Number of variables    = 11 + #parameters
  
  // @todo Please implement/augment if required

  double u = Q[1];
  double v = Q[2];
  double w = Q[3];


  double jacobian = Q[6];

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];

  
  
  F[0][ 0] = -jacobian*(q_x*u+q_y*v+q_z*w);
  F[0][ 1] = 0.0;
  F[0][ 2] = 0.0;
  F[0][ 3] = 0.0;

  F[1][ 0] = -jacobian*(r_x*u+r_y*v+r_z*w);
  F[1][ 1] = 0.0;
  F[1][ 2] = 0.0;
  F[1][ 3] = 0.0;

  F[2][ 0] = -jacobian*(s_x*u+s_y*v+s_z*w);
  F[2][ 1] = 0.0;
  F[2][ 2] = 0.0;
  F[2][ 3] = 0.0;
}


void Linear::MyLinearSolver::nonConservativeProduct(const double* const Q,const double* const gradQ,double* BgradQ) {

  double p_q = gradQ[0];  
  double u_q = gradQ[1];
  double v_q = gradQ[2];
  double w_q = gradQ[3];

  double p_r = gradQ[4];  
  double u_r = gradQ[5];
  double v_r = gradQ[6];
  double w_r = gradQ[7];

  double p_s = gradQ[8];  
  double u_s = gradQ[9];
  double v_s = gradQ[10];
  double w_s = gradQ[11];


  double jacobian = Q[6];

  double q_x=Q[7];
  double q_y=Q[8];
  double q_z=Q[9];

  double r_x=Q[10];
  double r_y=Q[11];
  double r_z=Q[12];

  double s_x=Q[13];
  double s_y=Q[14];
  double s_z=Q[15];
  
  BgradQ[0]= 0;  
  BgradQ[1]= -q_x*p_q;
  BgradQ[2]= -q_y*p_q;
  BgradQ[3]= -q_z*p_q;

  BgradQ[4]=0;  
  BgradQ[5]=-r_x*p_r;
  BgradQ[6]=-r_y*p_r;
  BgradQ[7]=-r_z*p_r;

  BgradQ[8]=0;  
  BgradQ[9] =-s_x*p_s;
  BgradQ[10]=-s_y*p_s;
  BgradQ[11]=-s_z*p_s;

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


void Linear::MyLinearSolver::pointSource(const double* const x,const double t,const double dt, double* forceVector, double* x0, int n){

  double pi = 3.14159265359;
  double sigma = 0.1149;
  double t0 = 0.7;
  double f = 0.0;
  double M0 = 1000.0;
  
  if(n == 0){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.25;
    x0[1] = 0.25;
    x0[2] = 0.25;
    
    forceVector[0] = 1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
  }else if(n == 1){
    f = M0*(1.0/(sigma*std::sqrt(2.0*pi)))*(std::exp(-((t-t0)*(t-t0))/(2.0*sigma*sigma)));
    
    x0[0] = 0.75;
    x0[1] = 0.75;
    x0[2] = 0.75;
    
    forceVector[0] = -1.*f;
    forceVector[1] = 0.0;
    forceVector[2] = 0.0;
    forceVector[3] = 0.0;
  }
    

}
void Linear::MyLinearSolver::multiplyMaterialParameterMatrix(const double* const Q, double* rhs){

  double rho = Q[4];  
  double c   = Q[5];
  double jacobian = Q[6];  
  double mu  = rho*c*c;

  rhs[0]=mu/jacobian * rhs[0];
  rhs[1]=1/rho * rhs[1];
  rhs[2]=1/rho * rhs[2];
  rhs[3]=1/rho * rhs[3];

  rhs[4]=mu/jacobian * rhs[4];
  rhs[5]=1/rho * rhs[5];
  rhs[6]=1/rho * rhs[6];
  rhs[7]=1/rho * rhs[7];

  rhs[8] =mu/jacobian * rhs[8];
  rhs[9] =1/rho * rhs[9];
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

  //std::cout<<isBoundaryFace<<std::endl;
  
  double norm_p_qr=1.0;
  double norm_m_qr=1.0;

  //std::cout<<isBoundaryFace<<std::endl;

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {

      double c_p=QR[idx_QLR(i,j,5)];
      double rho_p=QR[idx_QLR(i,j,4)];
      double mu_p=c_p*c_p*rho_p;

      double c_m=QL[idx_QLR(i,j,5)];
      double rho_m=QL[idx_QLR(i,j,4)];
      double mu_m=c_m*c_m*rho_m;

      double qm_x=QL[idx_QLR(i,j,7)];
      double qm_y=QL[idx_QLR(i,j,8)];
      double qm_z=QL[idx_QLR(i,j,9)];
      
      double rm_x=QL[idx_QLR(i,j,10)];
      double rm_y=QL[idx_QLR(i,j,11)];
      double rm_z=QL[idx_QLR(i,j,12)];
            
      double sm_x=QL[idx_QLR(i,j,13)];
      double sm_y=QL[idx_QLR(i,j,14)];
      double sm_z=QL[idx_QLR(i,j,15)];

      double qp_x=QR[idx_QLR(i,j,7)];
      double qp_y=QR[idx_QLR(i,j,8)];
      double qp_z=QR[idx_QLR(i,j,9)];
      
      double rp_x=QR[idx_QLR(i,j,10)];
      double rp_y=QR[idx_QLR(i,j,11)];
      double rp_z=QR[idx_QLR(i,j,12)];
            
      double sp_x=QR[idx_QLR(i,j,13)];
      double sp_y=QR[idx_QLR(i,j,14)];
      double sp_z=QR[idx_QLR(i,j,15)];


      
      if (normalNonZeroIndex == 0){
	
      	norm_m_qr = std::sqrt(qm_x*qm_x + qm_y*qm_y + qm_z*qm_z);
      	n_m[0] = qm_x/norm_m_qr;
      	n_m[1] = qm_y/norm_m_qr;
      	n_m[2] = qm_z/norm_m_qr;	

      	norm_p_qr = std::sqrt(qp_x*qp_x + qp_y*qp_y  + qp_z*qp_z );
      	n_p[0] = qp_x/norm_p_qr;
      	n_p[1] = qp_y/norm_p_qr;
      	n_p[2] = qp_z/norm_p_qr;	
      }
      
      if (normalNonZeroIndex == 1){

      	norm_m_qr = std::sqrt(rm_x*rm_x + rm_y*rm_y + rm_z*rm_z);
      	n_m[0] = rm_x/norm_m_qr;
      	n_m[1] = rm_y/norm_m_qr;
      	n_m[2] = rm_z/norm_m_qr;	

      	norm_p_qr = std::sqrt(rp_x*rp_x + rp_y*rp_y  + rp_z*rp_z);
      	n_p[0] = rp_x/norm_p_qr;
      	n_p[1] = rp_y/norm_p_qr;
      	n_p[2] = rp_z/norm_p_qr;	
      }

      if (normalNonZeroIndex == 2){

      	norm_m_qr = std::sqrt(sm_x*sm_x + sm_y*sm_y + sm_z*sm_z);
      	n_m[0] = sm_x/norm_m_qr;
      	n_m[1] = sm_y/norm_m_qr;
      	n_m[2] = sm_z/norm_m_qr;	

      	norm_p_qr = std::sqrt(sp_x*sp_x + sp_y*sp_y  + sp_z*sp_z);
      	n_p[0] = sp_x/norm_p_qr;
      	n_p[1] = sp_y/norm_p_qr;
      	n_p[2] = sp_z/norm_p_qr;	
      }


      // for(int k = 0 ; k< 3 ;k++){
      // 	n_m[k] = n[k];
      // 	n_p[k] = n[k];
      // }

      localBasis(n_p, m_p, l_p, 3);
      localBasis(n_m, m_m, l_m, 3);     
	
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
	
	  double r = 1.;
	
	  riemannSolver_BC0(v_p, sigma_p, z_p, r, v_hat_p, sigma_hat_p);
	  riemannSolver_BC0(v_m, sigma_m, z_m, r, v_hat_m, sigma_hat_m);
	
	}
      
	if (faceIndex == 3) {
	
	  double r = 0.;
	
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
