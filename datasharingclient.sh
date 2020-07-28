#!/bin/bash

cd code/

# Run protocol.
./bin/DataSharingClient 0 ../par/test.par.0.txt &
sleep 5
./bin/DataSharingClient 1 ../par/test.par.1.txt &
sleep 5
./bin/DataSharingClient 2 ../par/test.par.2.txt &
sleep 5
./bin/DataSharingClient 3 ../par/test.par.3.txt ../test_data/
wait