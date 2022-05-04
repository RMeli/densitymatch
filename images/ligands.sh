#!/bin/bash


for ligand in "BRD4" "CDK2"
do
    for sdf in $(ls ../experiments/ligands/${ligand}/*.sdf | grep -v "tran" | grep -v "murcko")
    do
        fin="${sdf%.*}_molgrid.pcd"
        fout="${sdf%.*}.pdf"
        fout=$(basename ${fout})

        python ../render_pcd.py ${fin} "ligands/${fout}"
    done
done
