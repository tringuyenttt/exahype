echo "Configure project for plenty-nodes scaling test (no output)."
mkdir plenty-nodes/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/plenty-nodes/Euler-no-output.exahype )
