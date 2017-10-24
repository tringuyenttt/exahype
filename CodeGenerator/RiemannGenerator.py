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
# Generates the code for the volume integral
# for a specific configuration
#


import Backend
import TemplatingUtils


class RiemannGenerator:
    m_context = {}
    
    # name of generated output file
    m_filename = "riemannSolver.cpp"


    def __init__(self, i_context):
        self.m_context = i_context


    def generateCode(self):        
        if(self.m_context['isLinear']):
        #TODO JMG implement correct linear case
            TemplatingUtils.renderAsFile('riemannSolverNonLinear_cpp.template', self.m_filename, self.m_context)
        else:            
            # render template
            TemplatingUtils.renderAsFile('riemannSolverNonLinear_cpp.template', self.m_filename, self.m_context)
