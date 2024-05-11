#!/bin/bash

###### USERS INPUT ############################################################

## Define the total time of the simulation and how often to output concentration field

TotalTime=0.1
runTimestep=1e-3
WriteTimestep=10

#### END OF USER INPUT #######################################################

cp system/controlDictTransport system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$WriteTimestep/g" system/controlDict
sed -i "s/runTimestep/$runTimestep/g" system/controlDict


if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    # Run multiSpeciesTransportFoam in parallel
    echo -e "Run multiSpeciesTransportFoam in parallel on $NP processors"
    mpirun -np $NP multiSpeciesTransportFoam -parallel  > multiSpeciesTransport.out
else
    #Run multiSpeciesTransportFoam
    echo -e "Run multiSpeciesTransportFoam"
    multiSpeciesTransportFoam > multiSpeciesTransport.out
fi
