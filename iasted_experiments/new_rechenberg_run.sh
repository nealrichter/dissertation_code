#!/bin/sh

mv -f data data.`date +%F-%T`
mkdir data

#onemax
./run_rechenberg_experiment.sh 1 1

#./run_rechenberg_experiment.sh 1 20

#./run_rechenberg_experiment.sh 1 100


#needle
./run_rechenberg_experiment.sh 2 1

#./run_rechenberg_experiment.sh 2 20

#./run_rechenberg_experiment.sh 2 100


#dectrap
./run_rechenberg_experiment.sh 3 1

#./run_rechenberg_experiment.sh 3 20

#./run_rechenberg_experiment.sh 3 100


#dec2trap
./run_rechenberg_experiment.sh 4 1

#./run_rechenberg_experiment.sh 4 20

#./run_rechenberg_experiment.sh 4 100


#2trap
./run_rechenberg_experiment.sh 5 1

#./run_rechenberg_experiment.sh 5 20

#./run_rechenberg_experiment.sh 5 100


#royal road
./run_rechenberg_experiment.sh 6 1

#./run_rechenberg_experiment.sh 6 20

#./run_rechenberg_experiment.sh 6 100

#long path
./run_rechenberg_experiment.sh 7 1

#./run_rechenberg_experiment.sh 7 20

#./run_rechenberg_experiment.sh 7 100

