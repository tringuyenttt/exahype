#!/bin/env python
##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# This file is pivotal to the code generator. It manages 
# the internal decisions about padding, triggers the code
# generation of the solver kernels and calls the back end 
# assembly code generation.
#

import os
import copy
import subprocess
import errno

import KernelsHeaderGenerator
import SpaceTimePredictorGenerator
import VolumeIntegralGenerator
import SurfaceIntegralGenerator
import RiemannGenerator
import SolutionUpdateGenerator
import AdjustSolutionGenerator
import StableTimeStepSizeGenerator
import WeightsGenerator
import DGMatrixGenerator
import ConfigurationParametersGenerator
import BoundaryConditionsGenerator
import ConverterGenerator
import GemmsCPPGenerator
import AMRRoutinesGenerator

g_config                 = {}
g_simdWidth              = {  'noarch' : 1,
                              'wsm'    : 2,
                              'snb'    : 4,
                              'hsw'    : 4,
                              'knc'    : 8,
                              'knl'    : 8 
                            }

def setConfig(i_config):
    global g_config
    g_config = i_config

def getSizeWithPadding(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = g_simdWidth[g_config['architecture']]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    return l_sizeWithPadding


def getPadWidth(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = g_simdWidth[g_config['architecture']]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    l_padWidth        = l_sizeWithPadding - i_sizeWithoutPadding
    return l_padWidth


def prepareOutputDirectory(i_outputDirectory):
    # create directory for output files if not existing
    try:
        os.makedirs(i_outputDirectory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    # remove all .cpp, .cpph, .c and .h files (we are in append mode!)
    for l_fileName in os.listdir(i_outputDirectory):
        _ , l_ext = os.path.splitext(l_fileName)
        if(l_ext in ['.cpp', '.cpph', '.c', '.h']):
            os.remove(i_outputDirectory + "/" + l_fileName)


def generateContext(i_config):
    context = copy.copy(i_config)
    context['nVarPad']  = getSizeWithPadding(context['nVar'])
    context['nParPad']  = getSizeWithPadding(context['nPar'])
    context['nDataPad'] = getSizeWithPadding(context['nData'])
    context['nDofPad']  = getSizeWithPadding(context['nDof'])
    context['nDof3D']   = 1 if context['nDim'] == 2 else context['nDof']
    context['isLinear'] = context['numerics'] == "linear"
    context['solverHeader']      = context['solverName'].split('::')[1] + '.h'
    context['codeNamespaceList'] = context['codeNamespace'].split('::')
    context['guardNamespace']    = '_'.join(context['codeNamespaceList']).upper()
    return context

    
def generateComputeKernels():
    kernelsHeaderGenerator = KernelsHeaderGenerator.KernelsHeaderGenerator(generateContext(g_config))
    kernelsHeaderGenerator.generateCode()
    spaceTimePredictorGenerator = SpaceTimePredictorGenerator.SpaceTimePredictorGenerator(generateContext(g_config))
    spaceTimePredictorGenerator.generateCode()
    volumeIntegralGenerator = VolumeIntegralGenerator.VolumeIntegralGenerator(generateContext(g_config))
    volumeIntegralGenerator.generateCode()
    surfaceIntegralGenerator = SurfaceIntegralGenerator.SurfaceIntegralGenerator(generateContext(g_config))
    surfaceIntegralGenerator.generateCode()
    riemannGenerator = RiemannGenerator.RiemannGenerator(generateContext(g_config))
    riemannGenerator.generateCode()
    solutionUpdateGenerator = SolutionUpdateGenerator.SolutionUpdateGenerator(generateContext(g_config))
    solutionUpdateGenerator.generateCode()
    adjustSolutionGenerator = AdjustSolutionGenerator.AdjustSolutionGenerator(generateContext(g_config))
    adjustSolutionGenerator.generateCode()
    stableTimeStepSizeGenerator = StableTimeStepSizeGenerator.StableTimeStepSizeGenerator(generateContext(g_config))
    stableTimeStepSizeGenerator.generateCode()
    weightsGenerator = WeightsGenerator.WeightsGenerator(generateContext(g_config))
    weightsGenerator.generateCode()
    dgMatrixGenerator = DGMatrixGenerator.DGMatrixGenerator(generateContext(g_config))
    dgMatrixGenerator.generateCode()
    gemmsCPPGenerator = GemmsCPPGenerator.GemmsCPPGenerator(generateContext(g_config))
    gemmsCPPGenerator.generateCode()
    configurationParametersGenerator = ConfigurationParametersGenerator.ConfigurationParametersGenerator(generateContext(g_config))
    configurationParametersGenerator.generateCode()
    boundaryConditionsGenerator = BoundaryConditionsGenerator.BoundaryConditionsGenerator(generateContext(g_config))
    boundaryConditionsGenerator.generateCode()
    converterGenerator = ConverterGenerator.ConverterGenerator(generateContext(g_config))
    converterGenerator.generateCode()
    amrRoutinesGenerator = AMRRoutinesGenerator.AMRRoutinesGenerator(generateContext(g_config))
    amrRoutinesGenerator.generateCode()


def executeLibxsmmGenerator(i_commandLineParameters):
    l_bashCommand = g_config['pathToLibxsmmGemmGenerator'] + i_commandLineParameters
    subprocess.call(l_bashCommand.split())


def generateAssemblerCode(i_outputFileName,
                          i_matmulConfigList):
    l_pathToAsmFile = os.path.splitext(i_outputFileName)[0]+'.c'
    for l_matmul in i_matmulConfigList:
        # for plain assembly code (rather than inline assembly) choose dense_asm
        l_commandLineArguments = ' ' + "dense"  + \
                                 ' ' + os.path.join(g_config['pathToOutputDirectory'],l_pathToAsmFile) + \
                                 ' ' + g_config['codeNamespace'] + '::' + l_matmul.baseroutinename + \
                                 ' ' + str(l_matmul.M) + \
                                 ' ' + str(l_matmul.N) + \
                                 ' ' + str(l_matmul.K) + \
                                 ' ' + str(l_matmul.LDA) + \
                                 ' ' + str(l_matmul.LDB) + \
                                 ' ' + str(l_matmul.LDC) + \
                                 ' ' + str(l_matmul.alpha) + \
                                 ' ' + str(l_matmul.beta) + \
                                 ' ' + str(l_matmul.alignment_A) + \
                                 ' ' + str(l_matmul.alignment_C) + \
                                 ' ' + g_config['architecture'] + \
                                 ' ' + l_matmul.prefetchStrategy+ \
                                 ' ' + 'DP' #always use double precision, 'SP' for single
        executeLibxsmmGenerator(l_commandLineArguments)
