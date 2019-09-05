#/bin/bash

subject=$1
data_dir=$2
region_list=$3
target_region=$4
atlas=$5
hemisphere=$6

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

while read region
do

    outDir=${data_dir}Connectopy/${atlas}/${subject}/
    outBase=${subject}.${hemisphere}.${region}.2.${target_region}.Evecs.func.gii
    simBase=${subject}.${hemisphere}.Eta2.${region}.2.${target_region}.mat
    outFile=${outDir}${outBase}

    if [ ! -f ${outFile} ]; then

        python ${eigenmaps} -s ${subject} \
        -l ${data_dir}Labels/${atlas}/${subject}.${hemisphere}.aparc.32k_fs_LR.label.gii \
        -sr ${region} \
        -sim ${outDir}${simBase} \
        -d ${outDir} \
        -o ${outBase} \
        -hemi ${hemisphere}

    fi
done <${region_list}
