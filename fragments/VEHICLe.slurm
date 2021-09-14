#!/bin/bash
#SBATCH --job-name=VEHICLe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --partition=short
#SBATCH --clusters=htc
#SBATCH --array=1001-2000
#SBATCH --output=slurm/%x.%A_%a.out
#SBATCH --time=00:05:00

DM="/data/biggin/lina3015/densitymatch/"
DEV="${DM}/development"

echo ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}

frag="VEHICLe/fragment_${SLURM_ARRAY_TASK_ID}.sdf"
outname=$(echo ${frag} | sed "s/.sdf/.pcd/g")

time singularity run --nv -B ${DM}:${DM} --app python ${DEV}/ligan.sif \
    ${DM}/molgrid_to_pcd.py ${frag} --ligmap ${DM}/files/ligmap -o ${outname}