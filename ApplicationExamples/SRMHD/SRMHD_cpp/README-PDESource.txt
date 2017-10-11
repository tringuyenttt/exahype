Status:

Started to implement MHD, but the source terms have to be stated somehow.

In the Dumbser/Zanotti Fortran framework, the function signature looks
like this:

SUBROUTINE PDESource(S,Q,x,time) 
  ! ...
  S = 0.
  S(9) = - DivCleaning_kappa * Q(9)
END SUBROUTINE PDESource

We should have exactly this signature also here in ExaHyPE.

Update: Since the coding week we now can do such stuff. However, for
the MHD we stick to DivCleaning_kappa=0.0 as MD does not like this term.
