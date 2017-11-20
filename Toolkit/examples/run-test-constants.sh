#!/bin/bash

TOOLKIT="java -jar ./Toolkit/dist/ExaHyPE.jar --not-interactive"

# code output directories
mkdir -p output/test-constants-aderdg

examples="Toolkit/examples"

cd  ../../

$TOOLKIT $examples/test-constants-aderdg.exahype
