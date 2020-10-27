#!/bin/sh
cd /afs/cern.ch/work/s/siluo/wprime/

ID=$(( $2 + 1 ))
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.20.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh
echo "Job SampleType: $1 - Sample Number: $ID has begun..."

# root -l -b "hypothesis.cc+($1,$ID,-2)"
# root -l -b "TemplateTest2d.cc+($1,$ID,0,-2)"
# root -l -b "TemplateTest2d.cc+($1,$ID,1,-2)"
# root -l -b "TemplateTest2d.cc+($1,$ID,10,-2)"
# root -l -b "TemplateTest2d.cc+($1,$ID,11,-2)"
# root -l -b "RunModule.cc+($1,$ID,0,-2)"
# root -l -b "LepTop.cc+($1,$ID,0,-2)"
root -l -b "MatchAll.cc+($1,$ID,0,-2)"

echo "Job SampleType: $1 - Sample Number: $ID is done..."
