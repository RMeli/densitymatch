#!/bin/bash

export root=${HOME}/Documents/git/densitymatch/
export ligroot=${root}/experiments/ligands/
export fragroot=${root}/fragments/
export sensaasroot=${HOME}/Documents/git/sensaas/
export PATH=${sensaasroot}:${PATH}

sensaas(){
    lig=$1

    mkdir -p ${lig}
    cd ${lig}

    # Run sensaas
    python ${sensaasroot}/meta-sensaas.py ${ligroot}/${lig}.sdf ${fragroot}/VEHICLe-good.sdf
}
export -f sensaas # Export function for GNU parallel

parallel -j 6 sensaas ::: "BRD4/ligand-1" "BRD4/ligand-2" "BRD4/ligand-3" "BRD4/ligand-4" "BRD4/ligand-5" "BRD4/ligand-6" "BRD4/ligand-7" "BRD4/ligand-8" "BRD4/ligand-9" "BRD4/ligand-10"
