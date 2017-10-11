#include "CurvilinearTransformation.h"

#include "kernels/KernelUtils.h"
#include "kernels/DGMatrices.h"


// void Linear::CurvilinearTransformation::test(){
//   return;
// }

void transFiniteInterpolation(int num_points_x, int num_points_y, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){


   double mesh_size_x = 1.0/(num_points_x-1);
   double mesh_size_y = 1.0/(num_points_y-1);

   double r;
   double q;
   kernels::idx2 id_xy(num_points_x,num_points_y);


  for(int j = 0 ; j < num_points_y ; j ++){
    printf("%f \n",top_bnd_y[j]);
} 
    printf("\n");
   for(int j =0 ; j < num_points_y ; j++) {
   for(int i =0 ; i < num_points_x ; i++) {

      	q = (i)*mesh_size_x;
        r = (j)*mesh_size_y;
        
        curvilinear_x[id_xy(i,j)] = (1-q)*left_bnd_x[j]+q*right_bnd_x[j]+(1-r)*bottom_bnd_x[i]+r*top_bnd_x[i]-
               (1-q)*(1-r)*left_bnd_x[0]-q*(1-r)*right_bnd_x[0]-r*(1-q)*top_bnd_x[0]-
               (r*q)*top_bnd_x[num_points_x-1];


        curvilinear_y[id_xy(i,j)] = (1-q)*left_bnd_y[j]+q*right_bnd_y[j]+(1-r)*bottom_bnd_y[i]+r*top_bnd_y[i]-
               (1-q)*(1-r)*left_bnd_y[0]-q*(1-r)*right_bnd_y[0]-r*(1-q)*top_bnd_y[0]-
               (r*q)*top_bnd_y[num_points_x-1];


   }
   }

  for(int j = 0 ; j < num_points_y ; j ++){
printf("%f \n",curvilinear_y[id_xy(j,num_points_y-1)]);
} 

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
      a_x=lagrangeBasis(x,orig_mesh_x,j,num_nodes);
      a_y=lagrangeBasis(y,orig_mesh_y,i,num_nodes);
      result += dest_mesh[id_xy(j,i)] * a_x*a_y;
    }
  }
}

void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){

  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
interpolate(kernels::gaussLegendreNodes[num_nodes-1][i],kernels::gaussLegendreNodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(i,j)]);
    }
  }
}

void computeDerivatives_x(int i, int j , double* values , int num_nodes, double& der_x){
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
        der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xy(n,j)];
      }
  //   }
  // }

}

void computeDerivatives_y (int i, int j , double* values , int num_nodes, double& der_y){
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
         der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xy(i,n)];
      }
  //   }
  // }
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
