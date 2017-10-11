Applications/SRHD
=================

This new application kernel was started by Vasco at Thu 23. Jun in order to merge the different
SRHD kernels around.

Caveats on the build system
---------------------------

You actually *can* rerun the toolkit on this, but make sure you recover the old Makefile:

    git checkout Makefile

If you are starting from a fresh checkout, notice that the build system currently
cannot correctly resolve dependencies, neither in C nor in Fortran. Thus, our single
remaining module Parameters.mod gets not created. Do this manually by invoking

   gfortran -c Parameters.f90

and running ./compile.sh afterwards.

A tool for FORTRAN to detect dependencies is https://github.com/ZedThree/fort_depend.py
A tool for C to detect dependencies is `makedepend`.