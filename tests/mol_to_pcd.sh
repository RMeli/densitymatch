#!/bin/bash

python transform.py ligands/BRD4/ligand-*.sdf
python transform.py ligands/CDK2/*.sdf

for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ligands/${ligand}/*_tran.sdf)
    do
        fout="${sdf%.*}.pcd"

        singularity run --nv --app python ../development/densitymatch.sif \
            ../molgrid_to_pcd.py ${sdf} -m ../files/ligmap -r 0.5 -o ${fout} -v
    done
done