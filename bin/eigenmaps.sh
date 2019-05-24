#/bin/bash

subject=$1
data_dir=$2
source_region=$3
target_region=$4
hemisphere=$5

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

python ${eigenmaps} -s ${subject} \
-l ${data_dir}Labels/Desikan/${subject}.${H}.aparc.32k_fs_LR.label.gii \
-sr ${source_region} \
-sim ${data_dir}Connectopy/Regional/${subject}/${source_region}/${subject}.${H}.Eta2.${source_region}.2.${target_region}.mat \
-d ${data_dir}Connectopy/Regional/${subject}/${source_region}/ \
-o ${subject}.${H}.${source_region}.2.${target_region}.Evecs.func.gii -hemi ${hemisphere}