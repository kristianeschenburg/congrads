#/bin/bash

subj=$1
region_file=$2
data_dir=$3
hemisphere=$4

CONGRAD_DIR=${GIT_DIR}/congrads
eta_script=${CONGRAD_DIR}/bin/eta2_regions.py

if [ ${hemisphere} == "L" ]; then 
    H="LEFT"
elif [ ${hemisphere} == "R" ]; then
    H="RIGHT"
else
    echo "Error: Incorrect hemisphere option." > logfile.log
    exit 125
fi

while read reg
do

    python ${eta_script} -s ${subj} -f ${data_dir}RestingState/${subj}.rfMRI_Z-Trans_merged_CORTEX_${H}.mat \
    -l ${data_dir}Labels/Desikan/${subj}.${hemisphere}.aparc.32k_fs_LR.label.gii \
    -sr ${reg} -d ${data_dir}Connectopy/${subj}/Regional/${reg} -bo ${reg}.2.brain -hemi ${hemisphere}

done < ${region_file}
