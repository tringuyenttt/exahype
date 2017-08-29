rm -r *.o cfiles.mk ffiles.mk kernels
echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/convergence/GRMHD_ADERDG.exahype )
mkdir -p convergence/results
