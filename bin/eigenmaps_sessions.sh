#/bin/bash

subject=$1
data_dir=$2
region_file=$3
atlas=$4
session=$5
ori=$6
hemisphere=$7

git_dir=/mnt/parcellator/parcellation/GitHub
CONGRAD_DIR=${git_dir}/congrads
eigenmaps=${CONGRAD_DIR}/bin/eigenmaps.py

if [ ${hemisphere} == "L" ]; then
    H="L"
elif [ ${hemisphere} == "R" ]; then
    H="R"
else
    echo "Error: Incorrect hemisphere option." > logfile.log
    exit 125
fi

while read source_region
do

    outDir=${data_dir}Connectopy/Regional/${subject}/Templated/${session}/${ori}/
    outBase=${subject}.${H}.${source_region}.2.brain.Evecs.func.gii
    simBase=${subject}.${H}.Eta2.${source_region}.2.brain.mat
    outFile=${outDir}${outBase}

    if [ ! -f ${outFile} ]; then

        python ${eigenmaps} -s ${subject} \
        -l ${data_dir}Labels/${atlas}/${hemisphere}.MaxProb.label.gii \
        -sr ${source_region} \
        -sim ${outDir}${simBase} \
        -d ${outDir} \
        -o ${outBase} \
        -hemi ${hemisphere}

    fi
done <${region_file}