echo "Configure project for multicore scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler/multicore/Euler-output.exahype )
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels