#!/bin/sh

mv -f data data.`date +%F-%T`
mkdir data


#dec2trap
./run_experiment.sh 4 1

#./run_experiment.sh 4 20

#./run_experiment.sh 4 100


#2trap
./run_experiment.sh 5 1

#./run_experiment.sh 5 20

#./run_experiment.sh 5 100

