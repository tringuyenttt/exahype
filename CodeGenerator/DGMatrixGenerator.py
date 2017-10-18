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


class DGMatrixGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "DGMatrices"
    
    # quadrature nodes and weights mapped onto [0,1]
    m_xGPN       = []
    m_wGPN       = []
    

    def __init__(self, i_context):
        self.m_context = i_context
        self.m_wGPN, self.m_xGPN = Utils.getGaussLegendre(self.m_context['nDof'])
        

    def generateCode(self):
        l_padSize = self.m_context['nDofPad'] - self.m_context['nDof']
        self.m_context['nDofPad_seq'] = range(self.m_context['nDofPad'])
        self.m_context['nDofPadTimesnDof_seq'] = range(self.m_context['nDofPad']*self.m_context['nDof'])

        # [FLCoeff 0...0]; [FRCoeff 0...0];
        # right now FLCoeff, FRCoeff no pad (gives no benefit w.r.t libxsmm)
        FLCoeff, _ = Utils.BaseFunc1d(0.0, self.m_xGPN, self.m_context['nDof']) #is also F0
        FRCoeff, _ = Utils.BaseFunc1d(1.0, self.m_xGPN, self.m_context['nDof'])
        self.m_context['FLCoeff'] = Utils.vectorPad(FLCoeff, l_padSize)
        self.m_context['FRCoeff'] = Utils.vectorPad(FRCoeff, l_padSize)
        
        # Matrices are stored in column major order (so the padding should be on the bottom rows)
        # [ Mat ]
        # [0...0]
        # [0...0]
        
        # Kxi
        Kxi = Utils.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_context['nDof'])
        self.m_context['Kxi'] = Utils.matrixPadAndFlatten_ColMajor(Kxi,l_padSize)

        # Kxi_T
        Kxi_T = Utils.matrixTranspose(Kxi)
        self.m_context['Kxi_T'] =Utils.matrixPadAndFlatten_ColMajor(Kxi_T,l_padSize)

        # iK1
        iK1 = Utils.matrixTranspose(Utils.matrixInverse(Utils.assembleK1(Kxi, self.m_xGPN, self.m_context['nDof'])))
        self.m_context['iK1'] = Utils.matrixPadAndFlatten_ColMajor(iK1,l_padSize)

        # dudx
        MM   = Utils.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_context['nDof'])
        dudx = Utils.assembleDiscreteDerivativeOperator(MM,Kxi)
        self.m_context['dudx'] = Utils.matrixPadAndFlatten_ColMajor(dudx,l_padSize)
        
        # dudx_T
        dudx_T = Utils.matrixTranspose(dudx)
        self.m_context['dudx_T'] = Utils.matrixPadAndFlatten_ColMajor(dudx_T,l_padSize)
        
        #fineGridProjector1d
        fineGridProjector1d_0 = Utils.assembleFineGridProjector1d(self.m_xGPN, 0, self.m_context['nDof'])
        fineGridProjector1d_1 = Utils.assembleFineGridProjector1d(self.m_xGPN, 1, self.m_context['nDof'])
        fineGridProjector1d_2 = Utils.assembleFineGridProjector1d(self.m_xGPN, 2, self.m_context['nDof'])
        self.m_context['fineGridProjector1d_0'] = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_0,l_padSize)
        self.m_context['fineGridProjector1d_1'] = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_1,l_padSize)
        self.m_context['fineGridProjector1d_2'] = Utils.matrixPadAndFlatten_ColMajor(fineGridProjector1d_2,l_padSize)
        self.m_context['fineGridProjector1d_T_0'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_0,l_padSize)
        self.m_context['fineGridProjector1d_T_1'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_1,l_padSize)
        self.m_context['fineGridProjector1d_T_2'] = Utils.matrixPadAndFlatten_RowMajor(fineGridProjector1d_2,l_padSize)
        
        #generate files 
        TemplatingUtils.renderAsFile('DGMatrices_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('DGMatrices_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)

        

        
    
    