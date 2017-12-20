echo "Configure project for plenty-nodes scaling test (no output)."
mkdir plenty-nodes/results
rm -r *.o cfiles.mk ffiles.mk cipofiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc2/Euler_ADERDG/plenty-nodes/Euler_ADERDG-no-output.exahype )
