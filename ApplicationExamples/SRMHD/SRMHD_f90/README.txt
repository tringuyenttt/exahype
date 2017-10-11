
The MHDf90 Application
======================

This application works like the MHD application and actually borrows 90% of the code.
However, it uses the Fortran ADER kernel instead of the C++ ADER kernel.

Thus it only can run for nDim=3 and polynomial order = 3.

It also MUST be compiled with MODE=Release , as soon as assertions are enabled, the code
fails.
