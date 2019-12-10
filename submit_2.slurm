#!/bin/bash
#SBATCH -J BG
#SBATCH -p stakeholder-4g
#SBATCH --time=72:00:00
#SBATCH -n1
#SBATCH --mem-per-cpu=4000
#SBATCH -a 1-175
#SBATCH -o slurmlog/Log_BG_%A-%a.out

echo "starting at `date` on `hostname`"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID" #ID in stats
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID" #ID in stats
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" #the index within in a job

# root -l -b "wprime.cpp+(isample,$SLURM_ARRAY_TASK_ID,-2)"
# root -l -b "hypothesis.cc+(2,$SLURM_ARRAY_TASK_ID,-2)"
# root -l -b "MatchTest.cc+(2,$SLURM_ARRAY_TASK_ID,-2)"
# root -l -b "MakePFile.cc+(2,$SLURM_ARRAY_TASK_ID,-2)"
root -l -b "GenMinimizeBestPermTest.cc+(2,$SLURM_ARRAY_TASK_ID,-2)"

echo "ended at `date` on `hostname`"
exit 0
