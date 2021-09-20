#!/bin/bash
#SBATCH --job-name=VEHICLe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --clusters=arc
##SBATCH --array=0-4999
##SBATCH --array=0-4867
#SBATCH --array=2933,4422
#SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=01:00:00

shift=0

DM="/data/biggin/lina3015/densitymatch/"
DEV="${DM}/development"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

# Shift SLURM_ARRAY_TARK_ID
# Needed to go around the limitations on the range of elements
SHIFTED=$((${SLURM_ARRAY_TASK_ID} + ${shift}))

echo ${SHIFTED}

frag="VEHICLe/fragment_${SHIFTED}.sdf"
outname=$(echo ${frag} | sed "s/.sdf/.pcd/g")

time singularity run -B ${DM}:${DM} --app python ${DEV}/ligan.sif \
    ${DM}/molgrid_to_pcd.py ${frag} --ligmap ${DM}/files/ligmap -o ${outname}