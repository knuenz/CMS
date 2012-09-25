#!/bin/sh

work=`pwd`

area=/afs/cern.ch/work/k/knuenz/CMSSW_4_2_3/src/ErnestOptimization

cd $area

eval `scramv1 runtime -sh`

rm -f stopjob$1

if [ ! -f variables$1.txt ]
    then
    echo variables$1.txt does not exist
    exit
fi

cp variables$1.txt $work/
cp optimizeSig.C $work/

cd $work

root -q -b -l optimizeSig.C\(\"$1\",\"$2\" >& optimizeSig${1}.log &

while((1))
  do
  sleep 30
  cp optimizeSig${1}.* $area/
  if [ -f $area/stopjob$1 ]
      then
      echo stopjob$1 found
      killall root.exe
      rm stopjob$1
      exit
  fi
done
