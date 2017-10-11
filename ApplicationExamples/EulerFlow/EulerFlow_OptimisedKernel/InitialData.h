#ifndef __INITIAL_DATA_EULERFLOW__
#define __INITIAL_DATA_EULERFLOW__

#include <cstdlib>
#include <iostream>
#include <string>

void InitialData(const double* const x, double* Q, double t);

void ShuVortex2D(const double* const x, double* Q, double t);

void MovingGauss2D(const double* const x, double* V, double t);

#endif /* __INITIAL_DATA_EULERFLOW__ */
