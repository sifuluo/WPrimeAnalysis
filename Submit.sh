#!/bin/sh
cd /afs/cern.ch/work/s/siluo/wprime/

ID=$(( $2 + 1 ))

echo "Job SampleType: $1 - Sample Number: $ID has begun..."

root -l -b "MakePFile.cc+($1,$ID,-2)"

echo "Job SampleType: $1 - Sample Number: $ID is done..."
