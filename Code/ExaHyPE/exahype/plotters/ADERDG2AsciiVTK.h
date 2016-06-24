/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_ASCII_VTK_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_ASCII_VTK_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

namespace exahype {
namespace plotters {
class ADERDG2AsciiVTK;
}
}

class exahype::plotters::ADERDG2AsciiVTK
    : public exahype::plotters::Plotter::Device {
 private:
  static int FileCounter;

  const std::string _filename;
  const int _order;
  const int _unknowns;

  tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*
      _patchWriter;
  tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter*
      _gridWriter;

  tarch::plotter::griddata::Writer::VertexDataWriter* _timeStampDataWriter;
  std::vector<tarch::plotter::griddata::Writer::VertexDataWriter*>
      _vertexDataWriter;

 public:
  ADERDG2AsciiVTK(const std::string& filename, int order, int unknowns);
  virtual ~ADERDG2AsciiVTK();

  virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);
};

#endif
