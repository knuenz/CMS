#!/bin/sh

#settings:
iterations=200
JOBNAME=AccEffCutUnRe70

make

cp GenScript GenScript${JOBNAME}

for scenario in 4;do

rap_=1
for pT_ in 6 7;do
./GenScript${JOBNAME} ${rap_}rap ${pT_}pT ${iterations}iter ${scenario}scen 
done

rap_=2
for pT_ in 6 7 8;do
./GenScript${JOBNAME} ${rap_}rap ${pT_}pT ${iterations}iter ${scenario}scen 
done

done

rm GenScript${JOBNAME}
