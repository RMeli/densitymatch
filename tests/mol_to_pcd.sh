#!/bin/bash


for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ligands/${ligand}/*.sdf)
    do
        fout="${sdf%.*}.pcd"

        singularity run --nv --app python ../containers/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -r 0.5 -o ${fout} -v
    done
done