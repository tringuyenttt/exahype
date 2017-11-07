#include "AbstractDummySolver_ADERDG.h"

constexpr int nVar = Dummy::AbstractDummySolver_ADERDG::NumberOfVariables;

inline void set(double* Q, double val) {
	for(int i=0; i<nVar; i++) Q[i] = val;
}

inline void initialdata(const double* const x, double* Q) {
    double x2 = (DIMENSIONS==3) ? x[2] : 0;
    double radius = std::sqrt( x[0]*x[0] + x[1]*x[1] + x2*x2 );
    
    set(Q, 0);
    Q[0] = radius;
}
