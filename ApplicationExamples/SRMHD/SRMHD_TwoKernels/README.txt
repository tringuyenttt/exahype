Neccessary Environment Variables
--------------------------------

Make sure you set

export MIXEDLANG=Yes
export MODE=Release

MIXEDLANG is needed so Fortran kernels get into ffiles.mk and are compiled.
Release is probably needed because Assertions do not work with any Fortran kernels for some reason.

