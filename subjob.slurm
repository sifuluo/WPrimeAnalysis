#!/bin/bash
#-n32
#SBATCH --mem-per-cpu=4000
#SBATCH -p stakeholder-4g

root -l -b wprime.cc+

echo "ended"
exit 0
