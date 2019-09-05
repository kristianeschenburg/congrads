#/bin/bash

subject=$1
data_dir=$2
region_list=$3
atlas=$4
hemisphere=$5

git_dir=/mnt/parcellator/parcellation/GitHub
CONGRAD_DIR=${git_dir}/congrads
eigenmaps=${CONGRAD_DIR}/bin/eigenmaps.py

if [ ${hemisphere} == "L" ]; then
    H="LEFT"
elif [ ${hemisphere} == "R" ]; then
    H="RIGHT"
else
    echo "Error: Incorrect hemisphere option." > logfile.log
    exit 125
fi

while read source_region
do 

    outDir=${data_dir}Connectopy/Templated/${atlas}/${subject}/
    outBase=${subject}.${hemisphere}.${source_region}.2.brain.Evecs.func.gii
    simBase=${subject}.${hemisphere}.Eta2.${source_region}.2.brain.mat
    outFile=${outDir}${outBase}

    if [ ! -f ${outFile} ]; then

        python ${eigenmaps} -s ${subject} \
        -l ${data_dir}Labels/${atlas}/${hemisphere}.100.aparc.32k_fs_LR.label.gii \
        -sr ${source_region} \
        -sim ${outDir}${simBase} \
        -d ${outDir} \
        -o ${outBase} \
        -hemi ${hemisphere}

    fi
done