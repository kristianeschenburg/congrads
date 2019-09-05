#/bin/bash

subj=$1
region_file=$2
data_dir=$3
atlas=$4
hemisphere=$5

git_dir=/mnt/parcellator/parcellation/GitHub
CONGRAD_DIR=${git_dir}/congrads
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

    out_dir=${data_dir}/Connectopy/Templated/${atlas}/${subj}/
    out_reg=${out_dir}${subj}.${hemisphere}.Eta2.${reg}.2.brain.mat

    if [ ! -f ${out_reg} ]; then
        echo ${reg}
        
	    python ${eta_script} -s ${subj} -f ${data_dir}/RestingState/${subj}.rfMRI_Z-Trans_merged.CORTEX_${H}.func.gii \
        -l ${data_dir}/Labels/${atlas}/${hemisphere}.100.aparc.32k_fs_LR.label.gii \
        -sr ${reg} -d ${out_dir} -bo ${reg}.2.brain -hemi ${hemisphere}
    fi

done < ${region_file}