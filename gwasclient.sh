#!/bin/bash

cd code/

# Run protocol.
./bin/GwasClient 0 ../par/test.par.0.txt &
sleep 5
./bin/GwasClient 1 ../par/test.par.1.txt &
sleep 5
./bin/GwasClient 2 ../par/test.par.2.txt &
wait