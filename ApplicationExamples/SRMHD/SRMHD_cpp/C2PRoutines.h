#ifndef __C2P_ROUTINES_CPP__
#define __C2P_ROUTINES_CPP__

/**
 * C definitions to use the MHD Con2Prim functions
 **/

extern "C" {

// FORTRAN functions called by C
void rtsafe_c2p_rmhd1_(double* X1,double* X2,double* XACC,double* gam,double* d,
    double* e,double* s2,double* b2,double* sb2,double* /* OUT */ w,bool* /* OUT */  failed);
void func_c2p_rmhd1_(double* x,double* /* OUT */ f,double* /* OUT */ df,double* gam,double* d,double* e,
    double* s2,double* b2,double* sb2,double* /* OUT */ w_out);


} /* extern "C" */
#endif /* __C2P_ROUTINES_CPP__ */
