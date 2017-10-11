# SVEC: State Vector Enhanced Code

SVEC is a code for defining the components of hyperbolic PDEs in a symbolic
manner. It defines and uses the tensish algebra-free tensor library and
some basic C++ templating.

The name emphazises the central nature of basically storing (linearized)
state vectors of the system and providing high level access.

Svec was written by SvenK in Sept 2017 for writing down the GRMHD equations
in the scope of ExaHyPE (exahype.eu).

## Concepts

- Stores and operates on a single spacetime-local state vector
- Strong attention on unneccessary copies: Shadowing (aliasing) of
  data structures instead of scatter/gather operations
- Use of a minimalistic algebra-free tensor library
- Readable PDEs with explicit contra- and covariant symbols
- Clean API towards ExaHyPE

## Upcoming Features

- Coupling to the Margherita microphysics code for advanced Cons2Prim
  operations and realistic equations of states.
- Coupling to the Einsteins Equations (FOCCZ4 PDE code)

## Copying and contact

Main authors:

  Sven Koeppel <koeppel@fias.uni-frankfurt.de>

The code is part of ExaHyPE and similiarly released under the
BSD 3 Open Source License. Copyright (c) 2017 http://exahype.eu,
All rights reserved.

The project has received funding from the European Union's Horizon
2020 research and innovation programme under grant agreement
No 671698. For copyrights and licensing, please consult the webpage.
