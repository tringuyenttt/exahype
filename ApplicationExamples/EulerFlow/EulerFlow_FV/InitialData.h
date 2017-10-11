#ifndef __InitialData_CLASS_HEADER__
#define __InitialData_CLASS_HEADER__

namespace Euler {

void rarefactionWave(const double* const x,double* Q);
void sodShockTube(const double* const x,double* Q);
void explosionProblem(const double* const x,double* Q);
void initialData(const double* const x,double* Q);

}

#endif // __InitialData_CLASS_HEADER__
