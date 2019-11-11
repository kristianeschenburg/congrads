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

if [ ${atlas} == 'Desikan' ]; then
	labext="aparc"
elif [ ${atlas} == 'Destrieux' ]; then
	labext="aparc.a2009s"
fi

while read source_region
do 

    outDir=${conn_dir}/${subject}/
    outBase=${subject}.${hemisphere}.${source_region}.2.brain.Evecs.Full
    simBase=${subject}.${hemisphere}.Eta2.${source_region}.2.brain.Full.mat
    outFile=${outDir}${outBase}

    # Full resting-state connectopy
    if [ ! -f ${outFile} ]; then

        python ${eigenmaps} -s ${subject} \
        -l ${labl_dir}/${atlas}/${hemisphere}.100.${labext}.32k_fs_LR.label.gii \
        -sr ${source_region} \
        -sim ${outDir}${simBase} \
        -o ${outFile} \
        -hemi ${hemisphere}

    fi

    # Loop over session connectopies
    for iter in $(seq 0 3); do
        
        outBase=${subject}.${hemisphere}.${source_region}.2.brain.Evecs.Iter.${iter}
        simBase=${subject}.${hemisphere}.Eta2.${source_region}.2.brain.Iter.${iter}.mat
        outFile=${outDir}${outBase}

        if [ ! -f ${outFile} ]; then

            python ${eigenmaps} -s ${subject} \
            -l ${labl_dir}/${atlas}/${hemisphere}.100.${labext}.32k_fs_LR.label.gii \
            -sr ${source_region} \
            -sim ${outDir}${simBase} \
            -o ${outFile} \
            -hemi ${hemisphere}

        fi
    done

done <${region_list}
