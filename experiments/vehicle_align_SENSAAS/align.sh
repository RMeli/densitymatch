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

parallel -j 4 sensaas ::: "CDK2/4fkp_B_LS5" "CDK2/4fkq_B_42K" "CDK2/4fkr_B_45K" "CDK2/4fks_B_46K" "CDK2/4fkt_B_48K" "CDK2/4fku_D_60K" "CDK2/4fkv_B_61K" "CDK2/4fkw_B_62K"

