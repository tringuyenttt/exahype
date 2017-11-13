#!/bin/bash
rename -v 's,tree-([0-9]+)-rank-([0-9]+).vtk,tree-rank-$2-$1.vtk,' tree*.vtk
