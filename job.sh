#!/bin/bash -l
#SBATCH --job-name=targets_make
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2-0
#SBATCH --mem=256G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.d.manurung@lumc.nl
#SBATCH --output=logs/targets.err

set -e
set -u
set -o pipefail

module purge
module load statistical/R/4.3.1

echo "Starting at `date`"

echo "Running on hosts: $SLURM_JOB_NODELIST"
echo "Running on $SLURM_JOB_NUM_NODES nodes."
echo "Running $SLURM_NTASKS tasks."
echo "Account: $SLURM_JOB_ACCOUNT"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node running script: $SLURMD_NODENAME"
echo "Submit host: $SLURM_SUBMIT_HOST"

R CMD BATCH run.R
rm -f .RData

echo "Program finished with exit code $? at: `date`"
