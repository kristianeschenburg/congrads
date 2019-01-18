#/bin/bash

subject=$1
source_region=$2
regions=$3

while read p
do

python eigenmaps.py -s ${subject} \
-r /mnt/parcellator/parcellation/HCP/Connectome_4/${subject}/Split_Surface_ROIS/Desikan_Killiany/func/L.${source_region}.func.gii \
-sim /mnt/parcellator/parcellation/parcellearning/Data/Connectopy/Regional/${subject}/${source_region}/${subject}.L.${source_region}.2.${p}.Eta2.mat \
-d /mnt/parcellator/parcellation/parcellearning/Data/Connectopy/Regional/${subject}/${source_region}/ \
-o ${subject}.L.${source_region}.2.${p}.Evecs.func.gii

done <${regions}