// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#ifndef __ComputeGlobalIntegralsLegendre_CLASS_HEADER__
#define __ComputeGlobalIntegralsLegendre_CLASS_HEADER__

#include "exahype/plotters/Plotter.h"
namespace Euler{
  class ComputeGlobalIntegralsLegendre;

  /**
   * Forward declaration
   */
  class MyEulerSolver;
}

/* I hope these modifications are not overwritten... */
#include "TimeSeriesReductions.h"
#include "MyEulerSolver.h"


#ifndef __NumberOfVariables__
static const int NumberOfVariables = Euler::MyEulerSolver::NumberOfVariables; // shortcut
#define __NumberOfVariables__
#endif

class Euler::ComputeGlobalIntegralsLegendre: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  private:
    reductions::MultipleReductionsWriter conserved;
    reductions::MultipleReductionsWriter primitives;
    reductions::MultipleReductionsWriter errors;
    reductions::ReductionsWriter statistics;
  public:
  ComputeGlobalIntegralsLegendre(MyEulerSolver&  solver);
  virtual ~ComputeGlobalIntegralsLegendre();
  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp) override;
};

#endif
