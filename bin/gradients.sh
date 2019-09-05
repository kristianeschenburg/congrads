#/bin/bash

subject=$1
conn_dir=$2
rois_dir=$3
surface_file=$4
region_file=$5
hemisphere=$6

echo ${region_file}

while read region
do

    metric=${conn_dir}${subject}.${hemisphere}.${region}.2.brain.Evecs.func.gii
    outGrad=${conn_dir}${subject}.${hemisphere}.${region}.2.brain.Gradient.func.gii
    outVecs=${conn_dir}${subject}.${hemisphere}.${region}.2.brain.GradientVectors.func.gii
    roi=${rois_dir}${hemisphere}.MaxProb.${region}.func.gii

    wb_command -metric-gradient ${surface_file} ${metric} ${outGrad}  \
    -vectors ${outVecs} \
    -roi ${roi}

done <${region_file}
