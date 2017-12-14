echo "Configure project for single-node scaling test (no output)."
mkdir single-node/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/coolmuc3_KNL/CCZ4/single-node/CCZ4-no-output.exahype )
