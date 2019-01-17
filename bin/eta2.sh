#/bin/bash

subject=$1
source_region=$2
regions=$3

while read p
do

python eta2.py -s ${subject} -f /mnt/parcellator/parcellation/parcellearning/Data/RestingState/${subject}.rfMRI_Z-Trans_merged_CORTEX_LEFT.mat -r /mnt/parcellator/parcellation/HCP/Connectome_4/${subject}/Split_Surface_ROIS/Desikan_Killiany/func/L.${source_region}.func.gii -t /mnt/parcellator/parcellation/HCP/Connectome_4/${subject}/Split_Surface_ROIS/Desikan_Killiany/func/L.${p}.func.gii -d /mnt/parcellator/parcellation/parcellearning/Data/Connectopy/Regional/${subject}/${source_region}/ -bo ${source_region}.2.${p}

done <${regions}
