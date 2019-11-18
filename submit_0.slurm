#!/bin/bash
#SBATCH -J FL
#SBATCH -p stakeholder-4g
#SBATCH --time=72:00:00
#SBATCH -n1
#SBATCH --mem-per-cpu=4000
#SBATCH -a 1-20
#SBATCH -o slurmlog/TDual_FormerLeptonic_log_%A-%a.out

echo "starting at `date` on `hostname`"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID" #ID in stats
echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID" #ID in stats
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" #the index within in a job

# root -l -b "wprime.cpp+(isample,$SLURM_ARRAY_TASK_ID)"
#root -l -b "hypothesis.cc+(0,$SLURM_ARRAY_TASK_ID)"
root -l -b "MatchTest.cc+(0,$SLURM_ARRAY_TASK_ID)"

echo "ended at `date` on `hostname`"
exit 0
