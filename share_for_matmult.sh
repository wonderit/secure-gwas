#!/bin/bash

cd code_ws/

# Run protocol.
./bin/DataSharingClient 0 ../par_ws/test.par.0.txt &
sleep 5
./bin/DataSharingClient 1 ../par_ws/test.par.1.txt &
sleep 5
./bin/DataSharingClient 2 ../par_ws/test.par.2.txt &
sleep 5
./bin/DataSharingClient 3 ../par_ws/test.par.3.txt ../test_data_ws_2/ &
wait