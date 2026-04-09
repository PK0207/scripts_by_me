#!/bin/bash
#SBATCH --job-name=line_fitting
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=128G
#SBATCH --account=ast7939
#SBATCH --qos=ast7939
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/fit_%j.out
#SBATCH --error=logs/fit_%j.err

#sinfo -o "%P %D %C %l %m" | sort -k5 -h  # sort by memory
#squeue -j <JOBID> -o "%R"
#sacctmgr show qos format=Name,GrpCPUs,GrpTRES,MaxCPUsPerUser -p | grep $(sacctmgr show user $USER -p | awk -F'|' '{print $2}' | tail -1)

module load python/3.10
source /home/pkottapalli/scripts_by_me/research24/astrovenv/bin/activate 

mkdir -p logs

CHUNK_SIZE=5

python -u sim_fitting.py $SLURM_CPUS_PER_TASK $CHUNK_SIZE
