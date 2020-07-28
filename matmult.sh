#!/bin/bash

cd code_ws/

# Run protocol.
./bin/MatMultClient 0 ../par_ws/test.par.0.txt &
sleep 5
./bin/MatMultClient 1 ../par_ws/test.par.1.txt &
sleep 5
./bin/MatMultClient 2 ../par_ws/test.par.2.txt &
wait