#/bin/bash

subj=$1
data_dir=$2
atlas=$3
hemisphere=$4

region_file=${data_dir}Labels/${atlas}/regions.txt

echo ${subj}
echo ${region_file}
echo ${data_dir}
echo ${atlas}
echo ${hemisphere}

git_dir=/mnt/parcellator/parcellation/GitHub
CONGRAD_DIR=${git_dir}/congrads
eta_script=${CONGRAD_DIR}/bin/eta2_regions.py

if [ ${atlas} == 'Desikan' ]; then
	labext='aparc'
elif [ ${atlas} == 'Destrieux' ]; then
	labext='aparc.a2009s'
fi

label_file=${data_dir}/Labels/${atlas}/${hemisphere}.100.${labext}.32k_fs_LR.label.gii
echo ${label_file}

feature_file=${data_dir}/RestingState/${subj}.${hemisphere}.rest.Z.merged.func.gii
echo ${feature_file}

out_dir=${data_dir}/Connectopy/Templated/${atlas}/${subj}/Pairwise/
mkdir -p ${out_dir}

while read sreg
do

    echo "Processing: "${sreg}

    out_sreg_dir=${out_dir}${sreg}/

    while read treg
    do

        out_reg=${out_sreg_dir}${subj}.${hemisphere}.Eta2.${sreg}.2.${treg}.Full.mat

        if [ ! -f ${out_sreg_dir} ]; then
            
            python ${eta_script} -s ${subj} -f ${feature_file} \
            -l ${label_file} \
            -sr ${sreg} -tr ${treg} -d ${out_sreg_dir} -bo ${sreg}.2.${treg} -hemi ${hemisphere} \
            -pi True -pf True
        fi

    done < ${region_file}
done < ${region_file}
