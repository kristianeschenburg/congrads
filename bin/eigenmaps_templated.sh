#/bin/bash

subject=$1
region_list=$2
atlas=$3
hemisphere=$4

data_dir=/mnt/parcellator/parcellation/parcellearning/Data/
conn_dir=${data_dir}Connectopy/Templated/${atlas}
labl_dir=${data_dir}Labels

git_dir=/mnt/parcellator/parcellation/GitHub
CONGRAD_DIR=${git_dir}/congrads
eigenmaps=${CONGRAD_DIR}/bin/eigenmaps.py

while read source_region
do 

    outDir=${conn_dir}/${subject}/
    outBase=${subject}.${hemisphere}.${source_region}.2.brain.Evecs
    simBase=${subject}.${hemisphere}.Eta2.${source_region}.2.brain.mat
    outFile=${outDir}${outBase}

    if [ ! -f ${outFile} ]; then

        python ${eigenmaps} -s ${subject} \
        -l ${labl_dir}/${atlas}/${hemisphere}.100.aparc.32k_fs_LR.label.gii \
        -sr ${source_region} \
        -sim ${outDir}${simBase} \
        -o ${outFile} \
        -hemi ${hemisphere}

    fi
done <${region_list}
