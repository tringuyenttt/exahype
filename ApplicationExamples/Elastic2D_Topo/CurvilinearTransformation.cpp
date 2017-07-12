#include "CurvilinearTransformation.h"

#include "kernels/KernelUtils.h"
#include "kernels/DGMatrices.h"


// void Linear::CurvilinearTransformation::test(){
//   return;
// }


void getBoundaryCurves(int num_points,double offset_x, double offset_y,double width_x, double width_y ,double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y){

 
  double pi = 3.14159265359;

  int nx = 1/width_x*num_points;
  int ny = 1/width_y*num_points;


  double dx= 1.0/(nx-1);
  double dy= 1.0/(ny-1);

  // int i_0 =  (offset_x-0.0)/width_x;
  // int j_0 =  (offset_y-0.0)/width_y;
  

  //std::cout<<i_0<<std::endl;
  //std::cout<<j_0<<std::endl;
  //std::cout<<std::endl;

  
  // for(int i = 0 ; i< num_points; i++){
  //   left_bnd_x[i] =  offset_x;
  //   right_bnd_x[i] = width_x+offset_x;
  //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
  //   top_bnd_x[i] = width_x*dx*i + offset_x;

  //   left_bnd_y[i] = width_y*dy*i + offset_y;
  //   right_bnd_y[i] = width_y*dy*i + offset_y;
  //   bottom_bnd_y[i] =offset_y;
  //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
  // }


  for(int i = 0 ; i< nx; i++){
    
    bottom_bnd_x[i] = 0 + dx*i;
    top_bnd_x[i]    = 0 + dx*i;

    
    bottom_bnd_y[i] = 0.0;
    top_bnd_y[i]    = 1 + 0.1*std::sin(2*pi*top_bnd_x[i]);

    //std::cout<<top_bnd_x[i]<<std::endl;
    //std::cout<<i<<std::endl;
    //std::cout<<std::endl;
  }

  //std::exit(-1);
  
  for(int i = 0 ; i< ny; i++){
    
    // left_bnd_x[i]  = 0;
    // right_bnd_x[i] = 1;
    
    left_bnd_x[i]  = 0;
    right_bnd_x[i] = 1;

    left_bnd_y[i]  = dy*i;
    right_bnd_y[i] = dy*i;
  }


  double h0y = (top_bnd_y[0] - bottom_bnd_y[0])/(ny-1);
  double hny = (top_bnd_y[nx-1] - bottom_bnd_y[nx-1])/(ny-1);

  for(int i = 0 ; i< ny; i++){
    left_bnd_y[i]  = bottom_bnd_y[0] + h0y*i;
    right_bnd_y[i] = bottom_bnd_y[nx-1] + hny*i;
  }
  
}

void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){


   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);

   double r;
   double q;
   kernels::idx2 id_xy(num_points,num_points);
   
   int j_0 = 0;
  
   for(int j =j_m ; j < j_p ; j++) {

     
     
     int i_0 = 0;
     
     for(int i =i_m ; i < i_p ; i++) {
       
       //std::cout<<i<<std::endl;
       
       q = (i)*mesh_size_x;
       r = (j)*mesh_size_y;
       
       curvilinear_x[id_xy(j_0,i_0)] = (1-q)*left_bnd_x[j]+q*right_bnd_x[j]+(1-r)*bottom_bnd_x[i]+r*top_bnd_x[i]-
	 (1-q)*(1-r)*left_bnd_x[0]-q*(1-r)*right_bnd_x[0]-r*(1-q)*top_bnd_x[0]-
	 (r*q)*top_bnd_x[mx-1];
       
       
       curvilinear_y[id_xy(j_0,i_0)] = (1-q)*left_bnd_y[j]+q*right_bnd_y[j]+(1-r)*bottom_bnd_y[i]+r*top_bnd_y[i]-
	 (1-q)*(1-r)*left_bnd_y[0]-q*(1-r)*right_bnd_y[0]-r*(1-q)*top_bnd_y[0]-
	 (r*q)*top_bnd_y[mx-1];
       
       i_0 = i_0+1; 
	 
       
       
     }
     //std::exit(-1);
     j_0 = j_0+1;
   }
   
   //  for(int j = 0 ; j < num_points_y ; j ++){
   // printf("%f \n",curvilinear_y[id_xy(j,num_points_y-1)]);
   //} 
   
}


double lagrangeBasis(double x,double* points,int i,int num_points){
  double result=1;
  for (int j = 0 ; j< num_points ; j ++){

    if (j != i) {
      result *= (x-points[j])/(points[i]-points[j]);
    }
  }
   
  return result;
}


void interpolate(double x, double y, double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes,double& result){

  double a_x=0;
  double a_y=0;
  
  result=0;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);

  
  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      a_x=lagrangeBasis(x,orig_mesh_x,i,num_nodes);
      a_y=lagrangeBasis(y,orig_mesh_y,j,num_nodes);
      result += dest_mesh[id_xy(j,i)] * a_x*a_y;
    }
  }
}

void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){

  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
      interpolate(kernels::gaussLegendreNodes[num_nodes-1][i],kernels::gaussLegendreNodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
    }
  }
}

void computeDerivatives_x(int i, int j , double* values , int num_nodes, double& der_x, double dx){
  //  double* vals_x;
  //  double* vals_y

    //  vals_x = malloc(sizeof(double)*num_nodes*num_nodes);
    //  vals_y = malloc(sizeof(double)*num_nodes*num_nodes);


  //  free(vals_x);
  //  free(vals_y);
  
   kernels::idx2 id_xy(num_nodes,num_nodes);

   der_x = 0.0;

  // for (int j = 0 ; j< num_nodes ; j ++){
  //   for (int i = 0 ; i< num_nodes ; i ++){
   for (int n = 0 ; n< num_nodes ; n ++){
     der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xy(j,n)]/dx;
   }
  //   }
  // }

}

void computeDerivatives_y (int i, int j , double* values , int num_nodes, double& der_y, double dy){
  //  double* vals_x;
  //  double* vals_y
    //  vals_x = malloc(sizeof(double)*num_nodes*num_nodes);
    //  vals_y = malloc(sizeof(double)*num_nodes*num_nodes);
  //  free(vals_x);
  //  free(vals_y);
  
  kernels::idx2 id_xy(num_nodes,num_nodes); 

  der_y = 0.0;
  // for (int j = 0 ; j< num_nodes ; j ++){
  //   for (int i = 0 ; i< num_nodes ; i ++){
      for (int n = 0 ; n< num_nodes ; n ++){
         der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xy(n,i)]/dy;
      }
  //   }
  // }
}



void metricDerivativesAndJacobian(int num_nodes, double* curvilinear_x, double* curvilinear_y,double* gl_vals_x,double* gl_vals_y,double* q_x,double* q_y,double* r_x,double* r_y,double* jacobian, double dx, double dy){

  double* unif_mesh= new double[num_nodes];
  for(int i = 0; i< num_nodes ; i++){
    unif_mesh[i]=i*1.0/(num_nodes-1);
  }

  getValuesAtQuadNodes(unif_mesh,unif_mesh,curvilinear_x,num_nodes,gl_vals_x);
  getValuesAtQuadNodes(unif_mesh,unif_mesh,curvilinear_y,num_nodes,gl_vals_y);

  double x_der_x; 
  double x_der_y;
  double y_der_x; 
  double y_der_y;
  
  kernels::idx2 id_xy(num_nodes,num_nodes);
  
  for(int j = 0 ; j < num_nodes ; j ++){
    for(int i = 0 ; i< num_nodes ; i ++){

      computeDerivatives_x(j,i,gl_vals_x,num_nodes,x_der_x, dx);
      computeDerivatives_y(j,i,gl_vals_x,num_nodes,x_der_y, dy);

      computeDerivatives_x(j,i,gl_vals_y,num_nodes,y_der_x, dx);
      computeDerivatives_y(j,i,gl_vals_y,num_nodes,y_der_y, dy);

// J = x_q.*y_r - y_q.*x_r; % Jacobian (determinant of metric)
// q_x = y_r./J;   % 1/J y_r
// r_x = - y_q./J; % -1/J y_q
// q_y = - x_r./J; % -1/J x_r
// r_y = x_q./J;   %  1/J x_q
		     
     jacobian[id_xy(j,i)]=x_der_x*y_der_y-x_der_y*y_der_x;
     
     q_x[id_xy(j,i)]=y_der_y/jacobian[id_xy(j,i)];
     q_y[id_xy(j,i)]=-x_der_y/jacobian[id_xy(j,i)];
     r_x[id_xy(j,i)]=-y_der_x/jacobian[id_xy(j,i)];
     r_y[id_xy(j,i)]=x_der_x/jacobian[id_xy(j,i)];

     //      printf("%f\n",der_x);
     //      printf("%f\n",der_y);
    }
 }

  


}


// void computeCurvilinearTransformationCoefficients(double* dest_mesh_x,double* dest_mesh_y, int num_nodes){

// double* dest_mesh_x_der_x;
// double* dest_mesh_x_der_y;
// double* dest_mesh_y_der_x;
// double* dest_mesh_y_der_y;

// dest_mesh_x_der_x= malloc(sizeof(double)*num_nodes*num_nodes);
// dest_mesh_x_der_y= malloc(sizeof(double)*num_nodes*num_nodes);
// dest_mesh_x_der_y= malloc(sizeof(double)*num_nodes*num_nodes);


// }
