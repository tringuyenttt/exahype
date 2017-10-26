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


import Backend
import TemplatingUtils
import Utils #matrix operation and build functions


class QuadratureGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "Quadrature"   

    def __init__(self, i_context):
        self.m_context = i_context
        
        

    def generateCode(self):
        if self.m_context['quadratureType'] == 'Gauss-Legendre':
            l_weights, l_nodes = Utils.getGaussLegendre(self.m_context['nDof'])
            
            #Nodes
            self.m_context['nodes_seq'] = range(self.m_context['nDof'])
            self.m_context['nodes'] = l_nodes
            
            #Weights
            l_weightsVector      = Utils.vectorPad(l_weights, self.m_context['nDofPad'] - self.m_context['nDof'])
            self.m_context['weights1'] = l_weightsVector
            self.m_context['w1Size'] = len(self.m_context['weights1'])
            self.m_context['w1_seq'] = range(self.m_context['w1Size'])
            if(self.m_context['nDim'] == 2):
                # weightsVector is wGPN itself
                l_weightsVector      = Utils.vectorPad(l_weights, Backend.getPadWidth(len(l_weights)))
                self.m_context['weights2'] = l_weightsVector
                self.m_context['w2Size'] = len(self.m_context['weights2'])
                self.m_context['w2_seq'] = range(self.m_context['w2Size'])

                # all combinations of two weights, written as an 1D array
                l_weightsVector = [l_weights[i] * l_weights[j] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof'])]
                l_weightsVector = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
                self.m_context['weights3'] = l_weightsVector
                self.m_context['w3Size'] = len(self.m_context['weights3'])
                self.m_context['w3_seq'] = range(self.m_context['w3Size'])

            elif(self.m_context['nDim'] == 3):
                # all combinations of two weights, written as an 1D array
                l_weightsVector = [l_weights[i] * l_weights[j] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof'])]
                l_weightsVector      = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
                self.m_context['weights2'] = l_weightsVector
                self.m_context['w2Size'] = len(self.m_context['weights2'])
                self.m_context['w2_seq'] = range(self.m_context['w2Size'])

                # all combination of three weights, written as an 1D array
                l_weightsVector = [l_weights[i] * l_weights[j] * l_weights[k] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof']) for k in range(self.m_context['nDof'])]
                l_weightsVector      = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
                self.m_context['weights3'] = l_weightsVector
                self.m_context['w3Size'] = len(self.m_context['weights3'])
                self.m_context['w3_seq'] = range(self.m_context['w3Size'])
                
            else:
                print("QuadratureGenerator.generateCode(): nDim == "+str(self.m_context['nDim'])+" not supported")

            
        else:
            print("QuadratureGenerator.generateCode(): quadratureType == "+str(self.m_context['quadratureType'])+" not supported")
            return #exit without generating 
        
        #generate files 
        TemplatingUtils.renderAsFile('Quadrature_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('Quadrature_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)
