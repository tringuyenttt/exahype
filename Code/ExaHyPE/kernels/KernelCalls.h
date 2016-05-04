// This file is part of the ExaHyPE project. For conditions of distribution and
// use, please see the copyright notice at www.exahype.eu

#ifndef _EXAHYPE_KERNELS_KERNEL_CALLS_H_
#define _EXAHYPE_KERNELS_KERNEL_CALLS_H_

#include "exahype/Parser.h"

namespace kernels {
// Is implemented within the application folder generated by the toolkit
void initSolvers(const exahype::Parser& parser);

/**
 * This callback enables the user to deallocate memory.
 */
void finalise();
}

#endif
