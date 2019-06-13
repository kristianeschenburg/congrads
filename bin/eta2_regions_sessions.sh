#/bin/bash

subj=$1
region_file=$2
data_dir=$3
atlas=$4
ori=$5
session=$6
hemisphere=$7

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

    python ${eta_script} -s ${subj} \
    -f ${data_dir}RestingState/Sessions/${subj}.rfMRI_REST${session}_${ori}_Z-Trans.CORTEX_${H}.func.gii \
    -l ${data_dir}Labels/${atlas}/${H}.MaxProb.label.gii \
    -sr ${reg} \
    -d ${data_dir}Connectopy/Regional/${subj}/Templated/${session}/${ori}/ \
    -bo ${reg}.2.brain \
    -hemi ${hemisphere}

done < ${region_file}
