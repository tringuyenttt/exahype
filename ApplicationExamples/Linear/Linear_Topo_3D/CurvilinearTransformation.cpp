#include "CurvilinearTransformation.h"

#include "kernels/KernelUtils.h"
#include "kernels/DGMatrices.h"
#if defined(_GLL)
#include "kernels/GaussLobattoQuadrature.h"
#else
#include "kernels/GaussLegendreQuadrature.h"
#endif


// void Linear::CurvilinearTransformation::test(){
//   return;
// }

double fault(double y, double z,double depth_y,double a_z, double b_z){
  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;
   std::cout << std::endl;  
   std::cout <<y <<std::endl;
   std::cout << 0.1*std::sin(2*pi*y) << std::endl;//*std::sin(2*pi*(z-(b_z+a_z)*0.5));
   //   return 0.1*std::sin(2*pi*y)*z*(1-z);//*std::sin(2*pi*(z-(b_z+a_z)*0.5));
   //std::cout <<(1.0/std::tan(angle)) * (y-depth_y/2.0) <<std::endl;
   //return (1.0/std::tan(angle)) * (y-depth_y/2.0);
   return z*(1-z);
}


double topography(double x, double z,double a_x, double b_x,double a_z, double b_z){

  double pi = 3.14159265359;
  double angle = 60.0/360.0 * 2 * pi;

  double Lx = b_x-a_x;
  double Lz = b_z-a_z;

  double topo;

  

  topo = (0.1*(x + z) + 3*(std::sin(4*pi*x/Lx+3.34)*std::cos(4*pi*x/Lx)
				* std::sin(4*pi*z/Lz+3.34)*std::cos(4*pi*z/Lz)));
  return topo;
}  


void getBoundaryCurves3D(int num_points,
			 double offset_x, double offset_y, double offset_z,
			 double width_x, double width_y , double width_z ,
			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
  double pi = 3.14159265359;

  int num_elt_x= std::ceil(1/width_x);
  int num_elt_y= std::ceil(1/width_y);
  int num_elt_z= std::ceil(1/width_z);

  int nx = num_elt_x*(num_points-1) + 1;
  int ny = num_elt_y*(num_points-1) + 1;
  int nz = num_elt_z*(num_points-1) + 1;  
  
  double dx= width_x/(num_points-1);
  double dy= width_y/(num_points-1);
  double dz= width_z/(num_points-1);

  // std::cout << nx  <<" " << ny  <<" "  << nz << std::endl;
  // std::cout << dx  <<" " << dy  <<" "  << dz << std::endl;

  // std::exit(-1);


  kernels::idx2 id_xy(ny,nx);// back front
  kernels::idx2 id_xz(nz,nx);// bottom top
  kernels::idx2 id_yz(nz,ny);//left right


  
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

  //top bottom
  for(int k = 0 ; k < nz ; k ++){
    for(int i = 0 ; i < nx ; i ++){
      //index = id_xz((elt_z*num_points)+k,(elt_x*num_points)+i);
      
      bottom_bnd_x[id_xz(k,i)] = dx*i;
      bottom_bnd_z[id_xz(k,i)] = dz*k;
      bottom_bnd_y[id_xz(k,i)] = 1;	    

      top_bnd_x[id_xz(k,i)] = dx*i;
      top_bnd_z[id_xz(k,i)] = dz*k;
      top_bnd_y[id_xz(k,i)] = 0+0.25*std::sin(2*pi*top_bnd_x[id_xz(k,i)])*std::sin(2*pi*top_bnd_z[id_xz(k,i)]);

	}
      }
    
  
  //left right
  for(int k = 0 ; k < nz ; k ++){
    for(int j = 0 ; j < ny ; j ++){
      //index = id_yz((elt_z*num_points)+k,(elt_y*num_points)+j);
      
      left_bnd_x[id_xz(k,j)]  = 0;	    
      left_bnd_y[id_xz(k,j)]  = dy*j;
      left_bnd_z[id_xz(k,j)]  = dz*k;
      
      right_bnd_x[id_xz(k,j)] = 1.0;
      right_bnd_y[id_xz(k,j)] = dy*j;
      right_bnd_z[id_xz(k,j)] = dz*k;
      
    }
  }
  
  //front back
  for(int j = 0 ; j < ny ; j ++){
    for(int i = 0 ; i < nx ; i ++){
      
      front_bnd_x[id_xy(j,i)]  = dx*i;
      front_bnd_y[id_xy(j,i)]  = dy*j;
      front_bnd_z[id_xy(j,i)]  = 0;
      
      back_bnd_x[id_xy(j,i)] = dx*i;
      back_bnd_y[id_xy(j,i)] = dy*j;
      back_bnd_z[id_xy(j,i)] = 1;
    }
  }
    
}

// void getBoundaryCurves3D(int num_points,
// 			 double offset_x, double offset_y, double offset_z,
// 			 double width_x, double width_y , double width_z ,
// 			 double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
// 			 double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
// 			 double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
// 			 double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
// 			 double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
// 			 double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
//   double pi = 3.14159265359;

//   int num_elt_x=1/width_x;
//   int num_elt_y=1/width_y;
//   int num_elt_z=1/width_z;

//   int nx = 1/width_x*num_points;
//   int ny = 1/width_y*num_points;
//   int nz = 1/width_z*num_points;  
  
//   double dx= width_x/(num_points-1);
//   double dy= width_y/(num_points-1);
//   double dz= width_z/(num_points-1);

//   std::cout << nx  <<" " << ny  <<" "  << nz << std::endl;
//   std::cout << dx  <<" " << dy  <<" "  << dz << std::endl;


//   kernels::idx2 id_xy(ny,nx);// back front
//   kernels::idx2 id_xz(nz,nx);// botton top
//   kernels::idx2 id_yz(nz,ny);//left right


  
//   // int i_0 =  (offset_x-0.0)/width_x;
//   // int j_0 =  (offset_y-0.0)/width_y;
  

//   //std::cout<<i_0<<std::endl;
//   //std::cout<<j_0<<std::endl;
//   //std::cout<<std::endl;

  
//   // for(int i = 0 ; i< num_points; i++){
//   //   left_bnd_x[i] =  offset_x;
//   //   right_bnd_x[i] = width_x+offset_x;
//   //   bottom_bnd_x[i] = width_x*dx*i + offset_x;
//   //   top_bnd_x[i] = width_x*dx*i + offset_x;

//   //   left_bnd_y[i] = width_y*dy*i + offset_y;
//   //   right_bnd_y[i] = width_y*dy*i + offset_y;
//   //   bottom_bnd_y[i] =offset_y;
//   //   top_bnd_y[i] = width_y+offset_y + 0.0*std::sin(2*pi*top_bnd_x[i]);
//   // }

//   //top bottom
//   int index ;
//   for(int elt_z = 0 ; elt_z< num_elt_z; elt_z++){
//     for(int elt_x = 0 ; elt_x< num_elt_x; elt_x++){
//       for(int k = 0 ; k < num_points ; k ++){
// 	for(int i = 0 ; i < num_points ; i ++){
// 	  index = id_xz((elt_z*num_points)+k,(elt_x*num_points)+i);

// 	  bottom_bnd_x[index] = width_x*elt_x+dx*i;
// 	  bottom_bnd_z[index] = width_z*elt_z+dz*k;
// 	  bottom_bnd_y[index] = 1;	    

// 	  top_bnd_x[index] = width_x*elt_x+dx*i;
// 	  top_bnd_z[index] = width_z*elt_z+dz*k;
// 	  top_bnd_y[index] = 0+0.0*std::sin(2*pi*top_bnd_x[index])*std::sin(2*pi*top_bnd_z[index]);

// 	  std::cout
// 	}
//       }
//     }
//   }
  
//   //left right
//   for(int elt_z = 0 ; elt_z< num_elt_z; elt_z++){
//     for(int elt_y = 0 ; elt_y< num_elt_y; elt_y++){
//       for(int k = 0 ; k < num_points ; k ++){
// 	for(int j = 0 ; j < num_points ; j ++){
// 	  index = id_yz((elt_z*num_points)+k,(elt_y*num_points)+j);

// 	  left_bnd_x[index]  = 0;	    
// 	  left_bnd_y[index]  = width_y*elt_y+dy*j;
// 	  left_bnd_z[index]  = width_z*elt_z+dz*k;
	  
// 	  right_bnd_x[index] = 1;
// 	  right_bnd_y[index] = width_y*elt_y+dy*j;
// 	  right_bnd_z[index] = width_z*elt_z+dz*k;
// 	}
//       }
//     }
//   }

//   //front back
//   for(int elt_y = 0 ; elt_y< num_elt_y; elt_y++){
//     for(int elt_x = 0 ; elt_x< num_elt_x; elt_x++){
//       for(int j = 0 ; j < num_points ; j ++){
// 	for(int i = 0 ; i < num_points ; i ++){
// 	  index = id_xy((elt_y*num_points)+j,(elt_x*num_points)+i);

// 	  front_bnd_x[index]  = width_x*elt_x+dx*i;
// 	  front_bnd_y[index]  = width_y*elt_y+dy*j;
// 	  front_bnd_z[index]  = 0;
	  
// 	  back_bnd_x[index] = width_x*elt_x+dx*i;
// 	  back_bnd_y[index] = width_y*elt_y+dy*j;
// 	  back_bnd_z[index] = 1;
// 	}
//       }
//     }
//   }


  
//   // double h0y = (top_bnd_y[0] - bottom_bnd_y[0])/(ny-1);
//   // double hny = (top_bnd_y[nx-1] - bottom_bnd_y[nx-1])/(ny-1);

//   // for(int i = 0 ; i< ny; i++){

//   //   left_bnd_y[i]  = bottom_bnd_y[0] + h0y*i;
//   //   right_bnd_y[i] = bottom_bnd_y[nx-1] + hny*i;
    
//   // }
  
// }


void getBoundaryCurves3D_fixedTopFace(int num_points,
				      double offset_x, double offset_y, double offset_z,
				      double width_x, double width_y , double width_z ,
				      double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
				      double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
				      double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
				      double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
				      double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
				      double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
  double pi = 3.14159265359;


  int nx = std::round(1/width_x)*(num_points-1) + 1;
  int ny = std::round(1/width_y)*(num_points-1) + 1;
  int nz = std::round(1/width_z)*(num_points-1) + 1;  

  double dx= 1.0/(nx-1);
  double dz= 1.0/(nz-1);
  double dy= 1.0/(ny-1);

  double depth=1.0;

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); //left right
  
  //(0.1*XT+1*sin(4*pi*XT/(xPMax-xPMin)+3.34).*cos(2*pi*(XT/(xPMax-xPMin)-0.5)+33.34))

  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){
      top_bnd_x[id_xz(k,i)] = 0+dx*i;      
      top_bnd_z[id_xz(k,i)] = 0+dz*k;
      top_bnd_y[id_xz(k,i)] = 0- 0*(0.1*(top_bnd_x[id_xz(k,i)] + top_bnd_z[id_xz(k,i)])
			       + 0.25*(std::sin(4*pi*top_bnd_x[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_x[id_xz(k,i)])
				       * std::sin(4*pi*top_bnd_z[id_xz(k,i)]+3.34)*std::cos(4*pi*top_bnd_z[id_xz(k,i)])));

    }
  }
  
  
  for(int i = 0 ; i< nx; i++){
    for(int k = 0 ; k< nz; k++){
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
      bottom_bnd_y[id_xz(k,i)] = depth;      
    }
  }


  getInterpolatedFace_fromBottomAndTop(nx,ny,nz,0,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       left_bnd_x,left_bnd_y,left_bnd_z);
    
  getInterpolatedFace_fromBottomAndTop(nx,ny,nz,1,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       right_bnd_x,right_bnd_y,right_bnd_z);
  
  getInterpolatedFace_fromBottomAndTop(nz,ny,nx,2,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       front_bnd_x,front_bnd_y,front_bnd_z);

  getInterpolatedFace_fromBottomAndTop(nz,ny,nx,3,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       back_bnd_x,back_bnd_y,back_bnd_z);
}



void getBoundaryCurves3D_fixedTopFace_forBlock(int num_points,
				      //				      double offset_x, double offset_y, double offset_z,
				      //				      double width_x, double width_y , double width_z ,
					       int nx, int ny, int nz, int n,
					       double width_x, double width_y , double width_z,	      
					       double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
					       double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
					       double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					       double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					       double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
					       double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

 
  
  double pi = 3.14159265359;


  double dx= 10*width_x/(nx-1);
  double dz= 10*width_y/(nz-1);
  double dy= 10*width_z/(ny-1);

  double depth=10.0;

  double x;
  double z;
  double y;

  kernels::idx2 id_xy(ny,nx); // back front
  kernels::idx2 id_xz(nz,nx); // botton top
  kernels::idx2 id_yz(nz,ny); // left right
  
  //(0.1*XT+1*sin(4*pi*XT/(xPMax-xPMin)+3.34).*cos(2*pi*(XT/(xPMax-xPMin)-0.5)+33.34))


  //given top surface
  for(int k = 0 ; k< nz; k++){  
    for(int i = 0 ; i< nx; i++){

      top_bnd_y[id_xz(k,i)] = 0;      
      top_bnd_z[id_xz(k,i)] = 0+dz*k;

     
      
      if(n == 0){
      top_bnd_x[id_xz(k,i)] = 0+dx*i;
      
      }else{
	top_bnd_x[id_xz(k,i)] = 5 + dx*i;
	
      }

      x = top_bnd_x[id_xz(k,i)];
      z = top_bnd_z[id_xz(k,i)];

      top_bnd_y[id_xz(k,i)] -= topography(x, z, 0.0, 10.0, 0.0, 10.0);
    }
  }
  
  for(int k = 0 ; k< nz; k++){
    for(int i = 0 ; i< nx; i++){
      bottom_bnd_y[id_xz(k,i)] = depth;      
      bottom_bnd_z[id_xz(k,i)] = top_bnd_z[id_xz(k,i)];
      bottom_bnd_x[id_xz(k,i)] = top_bnd_x[id_xz(k,i)];
      
    }
  }


  int n_left_right=ny;

  
  int n_block;
  int n_top_bottom;

  //left face
  n_block = nx;
  n_top_bottom = nz;

                                     //n_block, n_left_right, n_top_bottom,
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,0,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       left_bnd_x,left_bnd_y,left_bnd_z);

  //right face
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,1,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       right_bnd_x,right_bnd_y,right_bnd_z);


  //front face
  n_block = nz;
  n_top_bottom = nx;
  
  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,2,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       front_bnd_x,front_bnd_y,front_bnd_z);

  getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,3,
				       top_bnd_x,top_bnd_y,top_bnd_z,
				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
				       back_bnd_x,back_bnd_y,back_bnd_z);

  
}



// void getBoundaryCurves3D_cutOffTopography_withFault(int num_points,
// 				      //				      double offset_x, double offset_y, double offset_z,
// 				      //				      double width_x, double width_y , double width_z ,
// 						    int nx, int ny, int nz, int n,
// 						    double width_x, double width_y , double width_z,	      
// 						    double* left_bnd_x, double* left_bnd_y, double* left_bnd_z,
// 						    double* right_bnd_x, double* right_bnd_y, double* right_bnd_z,
// 						    double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
// 						    double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
// 						    double* front_bnd_x, double* front_bnd_y, double* front_bnd_z,
// 						    double* back_bnd_x, double* back_bnd_y, double* back_bnd_z){

  
//   double pi = 3.14159265359;


//   double dx= width_x/(nx-1);
//   double dz= width_y/(nz-1);
//   double dy= width_z/(ny-1);

//   double depth=1.0;

//   kernels::idx2 id_xy(ny,nx); // back front
//   kernels::idx2 id_xz(nz,nx); // botton top
//   kernels::idx2 id_yz(nz,ny); // left right




//   for(int k = 0 ; k< nz; k++){
//     for(jnt j = 0 ; j< ny; j++){
//       left_bnd_y[id_xy(k,j)] = 0.0 +dy*j
//       left_bnd_z[id_xy(k,j)] = 0.0 +dz*k

//       if(n==0){
// 	left_bnd_x[id_xy(k,j)] = 0.0;
//       }else{
// 	left_bnd_x[id_xy(k,j)] = 0.5-fault(left_bnd_y[id_xy(k,j)],left_bnd_z[id_xy(k,j)],1.0,0.0,1.0);
//       }
      
//     }
//   }


//   for(int k = 0 ; k< nz; k++){
//     for(jnt j = 0 ; j< ny; j++){
//       right_bnd_y[id_xy(k,j)] = 0.0 +dy*j
//       right_bnd_z[id_xy(k,j)] = 0.0 +dz*k

//       if(n==0){
// 	right_bnd_x[id_xy(k,j)] = 0.5-fault(right_bnd_y[id_xy(k,j)],right_bnd_z[id_xy(k,j)],1.0,0.0,1.0);
//       }else{
// 	right_bnd_x[id_xy(k,j)] = 1.0;
//       }
//     }
//   }


//   //top
//   double left_edge_x = new double[nx];
//   double right_edge_x = new double[nx];
//   double top_edge_x = new double[nx];
//   double bottom_edge_x = new double[nx];


//   double left_edge_z = new double[nz];
//   double right_edge_z = new double[nz];
//   double top_edge_z = new double[nz];
//   double bottom_edge_z = new double[nz];

//   // top face
//   for(int i_z=0; i_z < nz ;i_z++){
//     left_edge_x[i]  = left_bnd_x[id_yz(i_z,0)];
//     left_edge_z[i]  = left_bnd_z[id_yz(i_z,0)];

//     right_edge_x[i]  = right_bnd_x[id_yz(i_z,0)];
//     right_edge_z[i]  = right_bnd_z[id_yz(i_z,0)];    
//   }

//   double distance_x_bottom=(right_edge_x[0]-left_edge_x[0])/(nx-1);
//   double distance_x_top=(right_edge_x[nz-1]-left_edge_x[nz-1])/(nx-1);
  
//   for(int i_x=0; i_x < nx ;i_x++){
//     bottom_edge_x[i]  = left_edge_x[0] + distance_x_bottom*i_x;
//     bottom_edge_z[i]  = 0.0;

//     top_edge_x[i]  = left_edge_x[nz-1] + distance_x_top*i_x;
//     top_edge_z[i]  = 1.0;
//   }
  

//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_x,  right_edge_x,  top_edge_x,  bottom_edge_x,  top_bnd_x );
//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_z,  right_edge_z,  top_edge_z,  bottom_edge_z,  top_bnd_z );  
  
//   for(int k=0; k < nz ;k++){
//     for(int i=0; i < nx ;i++){
//       top_bnd_y[id_xz(k,i)]=0.0;
//     }
//   }


//   // bottom face
//   for(int i_z=0; i_z < nz ;i_z++){
//     left_edge_x[i]  = left_bnd_x[id_yz(i_z,ny-1)];
//     left_edge_z[i]  = left_bnd_z[id_yz(i_z,ny-1)];

//     right_edge_x[i]  = right_bnd_x[id_yz(i_z,ny-1)];
//     right_edge_z[i]  = right_bnd_z[id_yz(i_z,ny-1)];    
//   }

//   double distance_x_bottom=(right_edge_x[0]-left_edge_x[0])/(nx-1);
//   double distance_x_top=(right_edge_x[nz-1]-left_edge_x[nz-1])/(nx-1);
  
//   for(int i_x=0; i_x < nx ;i_x++){
//     bottom_edge_x[i]  = left_edge_x[0] + distance_x_bottom*i_x;
//     bottom_edge_z[i]  = 0.0;

//     top_edge_x[i]  = left_edge_x[nz-1] + distance_x_top*i_x;
//     top_edge_z[i]  = 1.0;
//   }
  

//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_x,  right_edge_x,  top_edge_x,  bottom_edge_x,  bottom_bnd_x );
//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_z,  right_edge_z,  top_edge_z,  bottom_edge_z,  bottom_bnd_z );  
  
//   for(int k=0; k < nz ;k++){
//     for(int i=0; i < nx ;i++){
//       bottom_bnd_y [id_xz(k,i)]=1.0;
//     }
//   }



//   // front face
//   for(int i_y=0; i_y < ny ;i_y++){
//     left_edge_x[i]  = left_bnd_x[id_yz(0,i_y)];
//     left_edge_y[i]  = left_bnd_y[id_yz(0,i_y)];

//     right_edge_x[i]  = right_bnd_x[id_yz(0,i_y)];
//     right_edge_z[i]  = right_bnd_z[id_yz(0,i_y)];    
//   }

//   double distance_x_bottom=(right_edge_x[0]-left_edge_x[0])/(nx-1);
//   double distance_x_top=(right_edge_x[nz-1]-left_edge_x[nz-1])/(nx-1);
  
//   for(int i_x=0; i_x < nx ;i_x++){
//     bottom_edge_x[i]  = left_edge_x[0] + distance_x_bottom*i_x;
//     bottom_edge_z[i]  = 0.0;

//     top_edge_x[i]  = left_edge_x[nz-1] + distance_x_top*i_x;
//     top_edge_z[i]  = 1.0;
//   }
  

//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_x,  right_edge_x,  top_edge_x,  bottom_edge_x,  front_bnd_x );
//   transFiniteInterpolation_singleCoordinate( nx, nz ,  left_edge_z,  right_edge_z,  top_edge_z,  bottom_edge_z,  front_bnd_z );  



  





//   for(int j = 0 ; j< ny; j++){
//     for(int i = 0 ; i< nx; i++){
//       front_bnd_y[id_xy(j,i)] = 0+dy*j;
//       front_bnd_z[id_xy(j,i)] = 0;

//       if(n == 0){
// 	front_bnd_x[id_xy(j,i)] = 0+dx*i;
//       }else{
// 	front_bnd_x[id_xy(j,i)] = 0.5+dx*i;
//       }
//     }
//   }

//   for(int j = 0 ; j< ny; j++){
//     for(int i = 0 ; i< nx; i++){
//       back_bnd_x[id_xy(j,i)] = front_bnd_x[id_xy(j,i)]
//       back_bnd_y[id_xy(j,i)] = front_bnd_y[id_xy(j,i)]
//       back_bnd_z[id_xy(j,i)] = 1.0;
//     }
//   }


//   int n_left_right=ny;

  
//   int n_block;
//   int n_top_bottom;

//   //left face
//   n_block = nx;
//   n_top_bottom = nz;

//                                      //n_block, n_left_right, n_top_bottom,
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,0,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       left_bnd_x,left_bnd_y,left_bnd_z);

//   //right face
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,1,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       right_bnd_x,right_bnd_y,right_bnd_z);


//   //front face
//   n_block = nz;
//   n_top_bottom = nx;
  
//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,2,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       front_bnd_x,front_bnd_y,front_bnd_z);

//   getInterpolatedFace_fromBottomAndTop(n_block, n_left_right, n_top_bottom,3,
// 				       top_bnd_x,top_bnd_y,top_bnd_z,
// 				       bottom_bnd_x,bottom_bnd_y,bottom_bnd_z,
// 				       back_bnd_x,back_bnd_y,back_bnd_z);


//   std::cout << "fault" << std::endl;
  
//   double z;
//   double y;
  
//   if(n==0){
//     for(int k = 0 ; k< nz; k++){
//       for(int i = 0 ; i< ny; i++){
// 	//	right_bnd_x[id_yz(k,i)] -= 1.0*fault(right_bnd_y[id_yz(k,i)],  right_bnd_z[id_yz(k,i)],depth,0,1);
// 	z=right_bnd_z[id_yz(i,k)];
// 	y=right_bnd_y[id_yz(i,k)];
// 	right_bnd_x[id_yz(k,i)] -=  z*(1-z)*(std::exp(-(y-0.5)*(y-0.5)/0.05));
//       }
//     }
//   }else{
//     for(int k = 0 ; k< nz; k++){
//       for(int i = 0 ; i< ny; i++){
// 	//	left_bnd_x[id_yz(k,i)] -= 1.0*fault(left_bnd_y[id_yz(k,i)],left_bnd_z[id_yz(k,i)],depth,0,1);
// 	z=left_bnd_z[id_yz(i,k)];
// 	y=left_bnd_y[id_yz(i,k)];
// 	left_bnd_x[id_yz(k,i)] -=   z*(1-z)*(std::exp(-(y-0.5)*(y-0.5)/0.05));

//       }
//     }
//   }
  
// }






void getInterpolatedFace_fromBottomAndTop( int n_block, int n_left_right, int n_top_bottom, int face,
					   double* top_bnd_x, double* top_bnd_y, double* top_bnd_z,
					   double* bottom_bnd_x, double* bottom_bnd_y, double* bottom_bnd_z,
					   double* face_x, double* face_y , double* face_z){

 
  double* top_edge_x= new double[n_top_bottom];
  double* top_edge_y= new double[n_top_bottom];
  double* top_edge_z= new double[n_top_bottom];

  double* bottom_edge_x= new double[n_top_bottom];
  double* bottom_edge_y= new double[n_top_bottom];
  double* bottom_edge_z= new double[n_top_bottom];

  double* left_edge_x= new double[n_left_right];
  double* left_edge_y= new double[n_left_right];
  double* left_edge_z= new double[n_left_right];

  double* right_edge_x= new double[n_left_right];
  double* right_edge_y= new double[n_left_right];
  double* right_edge_z= new double[n_left_right];

  // face 0:left 1:right 2:front 3:back
  bool left_right_face;
  bool first_face;
  switch(face){
  case(0): left_right_face=true; first_face=true ; break; //left
  case(1): left_right_face=true; first_face=false; break; //right
  case(2): left_right_face=false;first_face=true;  break; //bottom
  case(3): left_right_face=false;first_face=false; break; //top
  }


  int nx_on_top_bottom_face;
  int nz_on_top_bottom_face;
  if(left_right_face){
    nz_on_top_bottom_face=n_top_bottom;    
    nx_on_top_bottom_face=n_block;

  }else{
    nz_on_top_bottom_face=n_block;    
    nx_on_top_bottom_face=n_top_bottom;
  }
  
  kernels::idx2 id_xz(nz_on_top_bottom_face,nx_on_top_bottom_face);

  //generate top and bottom edge

  //  index = id_xz(nz-1-0,nx-1);
  // index = id_xz(0,nx-1);  
  // if(left_right & !first){
  //   std::cout <<  "top" << std::endl;
  //   std::cout <<  top_bnd_x[index] << std::endl;
  //   std::cout <<  top_bnd_y[index] << std::endl;;
  //   std::cout <<  top_bnd_z[index] << std::endl;;
  //   std::cout <<  "bottom" << std::endl;
  //   std::cout <<  bottom_bnd_x[index] << std::endl;
  //   std::cout <<  bottom_bnd_y[index] << std::endl;;
  //   std::cout <<  bottom_bnd_z[index] << std::endl;;     
  //   std::exit(-1);
  // } 


  // for(int i=0 ; i< nz; i++){
  //   for(int j=0 ; j< nx; j++){

  //     std::cout <<  "top" << std::endl;
  //     std::cout <<  top_bnd_x[id_xz(i,j)] << std::endl;
  //     std::cout <<  top_bnd_y[id_xz(i,j)] << std::endl;;
  //     std::cout <<  top_bnd_z[id_xz(i,j)] << std::endl;;
  //     std::cout <<  std::endl;
  //   }
  // }

  //   for(int i=0 ; i< nz; i++){
  //   for(int j=0 ; j< nx; j++){

  //     std::cout <<  "bottom" << std::endl;
  //     std::cout <<  bottom_bnd_x[id_xz(i,j)] << std::endl;
  //     std::cout <<  bottom_bnd_y[id_xz(i,j)] << std::endl;;
  //     std::cout <<  bottom_bnd_z[id_xz(i,j)] << std::endl;;
  //     std::cout <<  std::endl;
     
  //   }
  // } 
  
  //  std::cout << "top bottom edge" << std::endl;
  int index;
  for(int i=0; i < n_top_bottom ; i++){
    
    // left right front back
    if(left_right_face){
      if(first_face){ //left
	index = id_xz(i,0);
      }else{ //right
	index = id_xz(i,n_block-1);			
      }
    }else{//front or back
      if(first_face){ //front
	index=id_xz(0,i);
      }else{ //back
	//	index=id_xz(nz-1,nx-1-i);
	index=id_xz(n_block-1,i);
	// std::cout<< std::endl;	
	// std::cout << i << std::endl;
	// std::cout << index << std::endl;
      }
    }

    //    index = left_right ? (first ? id_xz(i,0) : id_xz(nz-1-i,nx-1)) : (first ? id_xz(0,i) :  id_xz(nz-1,nx-1-i));

    
    top_edge_x[i]=top_bnd_x[index];
    top_edge_y[i]=top_bnd_y[index];
    top_edge_z[i]=top_bnd_z[index];

    bottom_edge_x[i]=bottom_bnd_x[index];
    bottom_edge_y[i]=bottom_bnd_y[index];
    bottom_edge_z[i]=bottom_bnd_z[index];


    // std::cout << top_edge_x[i] <<std::endl;
    // std::cout << top_edge_y[i] <<std::endl;
    // std::cout << top_edge_z[i] <<std::endl;
                    
    // std::cout << bottom_edge_x[i] <<std::endl;
    // std::cout << bottom_edge_y[i] <<std::endl;
    // std::cout << bottom_edge_z[i] <<std::endl;
  }


  //
  double h_x=(bottom_edge_x[0]-top_edge_x[0])/(n_left_right-1);
  double h_y=(bottom_edge_y[0]-top_edge_y[0])/(n_left_right-1);
  double h_z=(bottom_edge_z[0]-top_edge_z[0])/(n_left_right-1);  

  

  // std::cout << "front edge" << std::endl;
  // std::cout << "hx: " << h_x<< "hy: " << h_y << "hz: " << h_z<<std::endl;
  for(int i=0; i < n_left_right ; i++){
    left_edge_x[i]=top_edge_x[0]+i*h_x;
    left_edge_y[i]=top_edge_y[0]+i*h_y;
    left_edge_z[i]=top_edge_z[0]+i*h_z;

    // std::cout << front_edge_x[i] <<std::endl;
    // std::cout << front_edge_y[i] <<std::endl;
    // std::cout << front_edge_z[i] <<std::endl;
  }

  h_x=(bottom_edge_x[n_top_bottom-1]-top_edge_x[n_top_bottom-1])/(n_left_right-1);
  h_y=(bottom_edge_y[n_top_bottom-1]-top_edge_y[n_top_bottom-1])/(n_left_right-1);
  h_z=(bottom_edge_z[n_top_bottom-1]-top_edge_z[n_top_bottom-1])/(n_left_right-1);  

  //  std::cout << "back edge" << std::endl;
  //  std::cout << "hx: " << h_x<< "hy: " << h_y << "hz: " << h_z<<std::endl;
  for(int i=0; i < n_left_right ; i++){
    right_edge_x[i]=top_edge_x[n_top_bottom-1]+i*h_x;
    right_edge_y[i]=top_edge_y[n_top_bottom-1]+i*h_y;
    right_edge_z[i]=top_edge_z[n_top_bottom-1]+i*h_z;

    // std::cout << back_edge_x[i] <<std::endl;
    // std::cout << back_edge_y[i] <<std::endl;
    // std::cout << back_edge_z[i] <<std::endl;

  }

  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_x,  right_edge_x,  top_edge_x,  bottom_edge_x,  face_x );
  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_y,  right_edge_y,  top_edge_y,  bottom_edge_y,  face_y );
  transFiniteInterpolation_singleCoordinate( n_top_bottom,  n_left_right,  left_edge_z,  right_edge_z,  top_edge_z,  bottom_edge_z,  face_z );


   // if(!left_right & !first){
   //  std::cout <<  "front" << std::endl;
   //  for(int i=0; i < ny ; i++){
   //  std::cout <<  front_edge_x[i] << std::endl;
   //  std::cout <<  front_edge_y[i] << std::endl;;
   //  std::cout <<  front_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;          
   //  }
   //  std::cout <<  "back" << std::endl;
   //  for(int i=0; i < ny ; i++){
   //    std::cout <<  back_edge_x[i] << std::endl;
   //    std::cout <<  back_edge_y[i] << std::endl;;
   //    std::cout <<  back_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;            
   //  }

   //  std::cout <<  "top" << std::endl;
   //  for(int i=0; i < nz ; i++){
   //    std::cout <<  top_edge_x[i] << std::endl;
   //    std::cout <<  top_edge_y[i] << std::endl;;
   //    std::cout <<  top_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;            
   //  }

   //  std::cout <<  "bottom" << std::endl;
   //  for(int i=0; i < nz ; i++){
   //    std::cout <<  bottom_edge_x[i] << std::endl;
   //    std::cout <<  bottom_edge_y[i] << std::endl;;
   //    std::cout <<  bottom_edge_z[i] << std::endl;;
   //    std::cout << std::endl;;      
   //  }

    

    //    std::exit(-1);
  //} 

  

  // kernels::idx2 id_zy(ny,nz);

  // std::cout << "face" <<std::endl;
  // std::cout << face <<std::endl;  


  // std::cout << "x" <<std::endl;    
  // for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_x[id_zy(j,i)] <<std::endl;
  //   }
  // }
   
  // std::cout << "y" <<std::endl;    
  //   for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_y[id_zy(j,i)] <<std::endl;
  //   }
  // }
    
  // std::cout << "z" <<std::endl;    
  //     for(int i=0; i < nz ; i++){
  //   for(int j=0; j < ny ; j++){
  //     std::cout << face_z[id_zy(j,i)] <<std::endl;
  //   }
  // }

}  


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
    top_bnd_y[i]    = 1 + 0.0*std::sin(2*pi*top_bnd_x[i]);

    //std::cout<<top_bnd_x[i]<<std::endl;
    //std::cout<<i<<std::endl;
    //std::cout<<std::endl;
  }

  //std::exit(-1);
  
  for(int i = 0 ; i< ny; i++){
    
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


void transFiniteInterpolation3D(int mx, int my, int mz,
				int k_m, int k_p ,
				int j_m, int j_p ,
				int i_m, int i_p ,
				int num_points,
				double width_x, double width_y, double width_z,
				double* left_bnd,
				double* right_bnd,
				double* bottom_bnd,
				double* top_bnd,
				double* front_bnd,
				double* back_bnd,
				double* curvilinear
				){
  kernels::idx3 id_xyz(num_points,num_points,num_points);
  kernels::idx2 id_xy(my,mx);// back front
  kernels::idx2 id_xz(mz,mx);// bottom top
  kernels::idx2 id_yz(mz,my);//left right

   // double mesh_size_x = width_x/(num_points-1);
   // double mesh_size_y = width_y/(num_points-1);
   // double mesh_size_z = width_z/(num_points-1);


  //double mesh_size_x = 1.0/(num_points-1);
  //double mesh_size_y = 1.0/(num_points-1);
  //double mesh_size_z = 1.0/(num_points-1);

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);
   double mesh_size_z = 1.0/(mz-1);

  // std::cout << mesh_size_x << " "<< mesh_size_y << " "<< mesh_size_z <<std::endl;
  // std::cout << mx << " "<< my << " "<< mz <<std::endl;

  int i_0;
  int j_0;
  int k_0;    

  double u,v,w;
  double uv,vw,uw,uvw1,uvw2;

  double q,r,s;

  for(int i=i_m ; i< i_p ; i++){
    i_0=i-i_m;
    q=mesh_size_x * i; //(i-std::floor((i+1)/num_points));
    //q=mesh_size_x *i_0;
    for(int j=j_m ; j< j_p ; j++){
      j_0=j-j_m;
      //r=mesh_size_y * (j-std::floor((j+1)/num_points));
      r=mesh_size_y *j;
      for(int k=k_m ; k< k_p ; k++){
	k_0=k-k_m;
	//s=mesh_size_z * (k-std::floor((k+1)/num_points));
	s=mesh_size_z * k;

	//std ::cout << i_0 << " "<< j_0 << " "<< k_0 <<std::endl;
	
	u=(1-q)*left_bnd[id_yz(k,j)]+q*right_bnd[id_yz(k,j)];

	v=(1-r)*top_bnd[id_xz(k,i)]+r*bottom_bnd[id_xz(k,i)];

	w=(1-s)*front_bnd[id_xy(j,i)]+s*back_bnd[id_xy(j,i)];

	uw=(1-q)* ((1 -r)*left_bnd[id_yz(k,0)] 
		   +r *left_bnd[id_yz(k,my-1)])
	  +q*((1-r)*right_bnd[id_yz(k,0)] 
	      +r *right_bnd[id_yz(k,my-1)]);

	uv=(1-r)*((1-s)*top_bnd[id_xz(0,i)]
		  + s*top_bnd[id_xz(mz-1,i)])
	  + r*((1-s)*bottom_bnd[id_xz(0,i)]
	       + s*bottom_bnd[id_xz(mz-1,i)]);

	vw = (1-s)*((1-q)*front_bnd[id_xy(j,0)] 
		    + q*front_bnd[id_xy(j,mx-1)])
	       + s*((1-q)*back_bnd[id_xy(j,0)] 
	              + q*back_bnd[id_xy(j,mx-1)]);
	
	uvw1=(1-q)*((1-r)*((1-s)*top_bnd[id_xz(0,0)] 
		              +s*top_bnd[id_xz(mz-1,0)])
		       +r*((1-s)*bottom_bnd[id_xz(0,0)]
			   +s*bottom_bnd[id_xz(mz-1,0)]));

	uvw2=q*((1-r)*((1-s)*top_bnd[id_xz(0,mx-1)] 
			 + s*top_bnd[id_xz(mz-1,mx-1)])
		+ r*((1-s)*bottom_bnd[id_xz(0,mx-1)]
		       + s*bottom_bnd[id_xz(mz-1,mx-1)]));

	curvilinear[id_xyz(k_0,j_0,i_0)] = u + v + w - uv - uw - vw + uvw1 + uvw2;


	// 	uw=(1-q)* ((1 -r)*left_bnd[id_yz(k,j_m)] 
	// 	   +r *left_bnd[id_yz(k,j_p-1)])
	//   +q*((1-r)*right_bnd[id_yz(k,j_m)] 
	//       +r *right_bnd[id_yz(k,j_p-1)]);

	// uv=(1-r)*((1-s)*top_bnd[id_xz(k_m,i)]
	// 	  + s*top_bnd[id_xz(k_p-1,i)])
	//   + r*((1-s)*bottom_bnd[id_xz(k_m,i)]
	//        + s*bottom_bnd[id_xz(k_p-1,i)]);

	// vw = (1-s)*((1-q)*front_bnd[id_xy(j,i_m)] 
	// 	    + q*front_bnd[id_xy(j,i_p-1)])
	//        + s*((1-q)*back_bnd[id_xy(j,i_m)] 
	//               + q*back_bnd[id_xy(j,i_p-1)]);
	
	// uvw1=(1-q)*((1-r)*((1-s)*top_bnd[id_xz(k_m,i_m)] 
	// 	              +s*top_bnd[id_xz(k_p-1,i_m)])
	// 	       +r*((1-s)*bottom_bnd[id_xz(k_m,i_m)]
	// 		   +s*bottom_bnd[id_xz(k_p-1,i_m)]));

	// uvw2=q*((1-r)*((1-s)*top_bnd[id_xz(k_m,i_p-1)] 
	// 		 + s*top_bnd[id_xz(k_p-1,i_p-1)])
	// 	+ r*((1-s)*bottom_bnd[id_xz(k_m,i_p-1)]
	// 	       + s*bottom_bnd[id_xz(k_p-1,i_p-1)]));

	// curvilinear[id_xyz(k_0,j_0,i_0)] = u + v + w - uv - uw - vw + uvw1 + uvw2;

      }
    }
  }
}


void transFiniteInterpolation_singleCoordinate(int mx, int my, double* left_bnd, double* right_bnd, double* bottom_bnd, double* top_bnd, double* curvilinear ){

   double mesh_size_x = 1.0/(mx-1);
   double mesh_size_y = 1.0/(my-1);
   
   // double mesh_size_x = width_x/(num_points-1);
   // double mesh_size_y = width_y/(num_points-1);
   // double mesh_size_z = width_z/(num_points-1);
  

   double r;
   double q;

   
   // std::cout << mesh_size_x << std::endl;
   // std::cout << mesh_size_y << std::endl;

   // std::cout << my << std::endl;
   // std::cout << mx << std::endl;

   
   //   std::exit(-1);

   kernels::idx2 id_xy(my,mx);
   
   for(int j =0 ; j < my ; j++) {
     for(int i = 0 ; i < mx ; i++) {
       
       q = (i)*mesh_size_x;
       r = (j)*mesh_size_y;
       
       curvilinear[id_xy(j,i)] = (1-q)*left_bnd[j]+q*right_bnd[j]+(1-r)*bottom_bnd[i]+r*top_bnd[i]-
	 (1-q)*(1-r)*left_bnd[0]-q*(1-r)*right_bnd[0]-r*(1-q)*top_bnd[0]-
	 (r*q)*top_bnd[mx-1];
     }
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
     j_0 = j_0+1;
   }
}

// void transFiniteInterpolation(int mx, int my, int j_m, int j_p, int i_m, int i_p, int num_points, double* left_bnd_x, double* left_bnd_y, double* right_bnd_x, double* right_bnd_y, double* bottom_bnd_x, double* bottom_bnd_y, double* top_bnd_x, double* top_bnd_y, double* curvilinear_x , double* curvilinear_y ){
  
//   transFiniteInterpolation_singleCoordinate(mx, my, left_bnd_x,right_bnd_x,bottom_bnd_x,top_bnd_x,curvilinear_x );
//   transFiniteInterpolation_singleCoordinate(mx, my, left_bnd_y,right_bnd_y,bottom_bnd_y,top_bnd_y,curvilinear_y );
// }


double lagrangeBasis(double x,double* points,int i,int num_points){
  double result=1;
  for (int j = 0 ; j< num_points ; j ++){

    if (j != i) {
      result *= (x-points[j])/(points[i]-points[j]);
    }
  }
   
  return result;
}



void interpolate3D(double x, double y, double z, double* orig_mesh_x , double* orig_mesh_y, double* orig_mesh_z, double* dest_mesh, int num_nodes,double& result){

  double a_x=0;
  double a_y=0;
  double a_z=0;  
  
  result=0;
  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);



  for (int k = 0 ; k< num_nodes ; k++){    
    for (int j = 0 ; j< num_nodes ; j ++){
      for (int i = 0 ; i< num_nodes ; i ++){
	a_x=lagrangeBasis(x,orig_mesh_x,i,num_nodes);
	a_y=lagrangeBasis(y,orig_mesh_y,j,num_nodes);
	a_z=lagrangeBasis(z,orig_mesh_z,k,num_nodes);
	result += dest_mesh[id_xyz(k,j,i)] * a_x*a_y*a_z;
      }
    }
  }
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



void getValuesAtQuadNodes3D(double* orig_mesh_x , double* orig_mesh_y, double* orig_mesh_z, double* dest_mesh, int num_nodes, double* results){

  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);


  for (int k = 0 ; k< num_nodes ; k ++){
    for (int j = 0 ; j< num_nodes ; j ++){
      for (int i = 0 ; i< num_nodes ; i ++){
#if defined(_GLL)
	interpolate3D(kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-i],kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-j],kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-k],orig_mesh_x,orig_mesh_y,orig_mesh_z,dest_mesh,num_nodes,results[id_xyz(k,j,i)]);
	#else
	interpolate3D(kernels::gaussLegendreNodes[num_nodes-1][i],kernels::gaussLegendreNodes[num_nodes-1][j],kernels::gaussLegendreNodes[num_nodes-1][k],orig_mesh_x,orig_mesh_y,orig_mesh_z,dest_mesh,num_nodes,results[id_xyz(k,j,i)]);
	#endif
      }
    }
  }
}



void getValuesAtQuadNodes(double* orig_mesh_x , double* orig_mesh_y, double* dest_mesh, int num_nodes, double* results){

  kernels::idx2 id_xy(num_nodes,num_nodes);

  for (int j = 0 ; j< num_nodes ; j ++){
    for (int i = 0 ; i< num_nodes ; i ++){
#if defined(_GLL)
      interpolate(kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-i],kernels::gaussLobattoNodes[num_nodes-1][num_nodes-1-j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
      #else
      interpolate(kernels::gaussLegendreNodes[num_nodes-1][i],kernels::gaussLegendreNodes[num_nodes-1][j],orig_mesh_x,orig_mesh_y,dest_mesh,num_nodes,results[id_xy(j,i)]);
#endif
    }
  }
}

void computeDerivatives_x(int i, int j , double* values , int num_nodes, double& der_x, double dx){
   kernels::idx2 id_xy(num_nodes,num_nodes);

   der_x = 0.0;

   for (int n = 0 ; n< num_nodes ; n ++){
     der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xy(j,n)]/dx;
   }
}

void computeDerivatives_y (int i, int j , double* values , int num_nodes, double& der_y, double dy){
  kernels::idx2 id_xy(num_nodes,num_nodes); 

  der_y = 0.0;
  for (int n = 0 ; n< num_nodes ; n ++){
    der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xy(n,i)]/dy;
  }
}


void computeDerivatives_x_3D(int i, int j , int k, double* values , int num_nodes, double& der_x, double dx){
  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);

  der_x = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_x += kernels::dudx[num_nodes-1][i][n] * values[id_xyz(k,j,n)]/dx;
  }

}

void computeDerivatives_y_3D (int i, int j,int k, double* values , int num_nodes, double& der_y, double dy){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_y = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_y += kernels::dudx[num_nodes-1][j][n] * values[id_xyz(k,n,i)]/dy;
  }
}

void computeDerivatives_z_3D (int i, int j , int k ,double* values , int num_nodes, double& der_z, double dz){
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes); 

  der_z = 0.0;

  for (int n = 0 ; n< num_nodes ; n ++){
    der_z += kernels::dudx[num_nodes-1][k][n] * values[id_xyz(n,j,i)]/dz;
  }

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


void metricDerivativesAndJacobian3D(int num_nodes,
				  double* curvilinear_x, double* curvilinear_y, double* curvilinear_z,
				  double* gl_vals_x, double* gl_vals_y, double* gl_vals_z,
				  double* q_x, double* q_y, double* q_z,
				  double* r_x, double* r_y, double* r_z,
				  double* s_x, double* s_y, double* s_z,				  
				  double* jacobian,
				  double dx, double dy, double dz
				  ){

  double* unif_mesh= new double[num_nodes];

  for(int i = 0; i< num_nodes ; i++){
    unif_mesh[i]=i*1.0/(num_nodes-1);
  }

  getValuesAtQuadNodes3D(unif_mesh,unif_mesh,unif_mesh,curvilinear_x,num_nodes,gl_vals_x);
  getValuesAtQuadNodes3D(unif_mesh,unif_mesh,unif_mesh,curvilinear_y,num_nodes,gl_vals_y);
  getValuesAtQuadNodes3D(unif_mesh,unif_mesh,unif_mesh,curvilinear_z,num_nodes,gl_vals_z);  

  double x_der_x; 
  double x_der_y;
  double x_der_z;  
  
  double y_der_x; 
  double y_der_y;
  double y_der_z;

  double z_der_x; 
  double z_der_y;
  double z_der_z;  

  
  kernels::idx3 id_xyz(num_nodes,num_nodes,num_nodes);
  
  for(int k = 0 ; k < num_nodes ; k ++){
    for(int j = 0 ; j < num_nodes ; j ++){
      for(int i = 0 ; i< num_nodes ; i ++){

  	computeDerivatives_x_3D(i,j,k,gl_vals_x,num_nodes,x_der_x, dx);
  	computeDerivatives_y_3D(i,j,k,gl_vals_x,num_nodes,x_der_y, dy);
  	computeDerivatives_z_3D(i,j,k,gl_vals_x,num_nodes,x_der_z, dz);

	computeDerivatives_x_3D(i,j,k,gl_vals_y,num_nodes,y_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_y,num_nodes,y_der_y, dy);
	computeDerivatives_z_3D(i,j,k,gl_vals_y,num_nodes,y_der_z, dz);

	computeDerivatives_x_3D(i,j,k,gl_vals_z,num_nodes,z_der_x, dx);
	computeDerivatives_y_3D(i,j,k,gl_vals_z,num_nodes,z_der_y, dy);
	computeDerivatives_z_3D(i,j,k,gl_vals_z,num_nodes,z_der_z, dz);
	
	jacobian[id_xyz(k,j,i)]=x_der_x*(y_der_y*z_der_z-y_der_z*z_der_y)
	  -x_der_y*(y_der_x*z_der_z-y_der_z*z_der_x)
	  +x_der_z*(y_der_x*z_der_y-y_der_y*z_der_x);
	
	q_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_y*z_der_z - z_der_y*y_der_z);
	r_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_x*y_der_z - y_der_x*z_der_z);
	s_x[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_x*z_der_y - z_der_x*y_der_y);
	
	q_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_y*x_der_z - x_der_y*z_der_z);
	r_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_x*z_der_z - z_der_x*x_der_z);
	s_y[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(z_der_x*x_der_y - x_der_x*z_der_y);
	
	q_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_y*y_der_z - y_der_y*x_der_z);
	r_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(y_der_x*x_der_z - x_der_x*y_der_z);
	s_z[id_xyz(k,j,i)] = (1.0/jacobian[id_xyz(k,j,i)])*(x_der_x*y_der_y - y_der_x*x_der_y);

      }
    }
  }

  


}

