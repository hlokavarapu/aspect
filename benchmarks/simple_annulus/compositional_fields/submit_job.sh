#!/bin/bash

interpolation=('bilinear' 'cell_average')
grid_resolution=('2' '3' '4' '5' '6')
job_template_filename='job_template.sh'

top_level_dir=`pwd`

for interpolation_scheme in ${interpolation[@]}; do
  cd $interpolation_scheme
  for grid_res in ${grid_resolution[@]}; do
    cd $grid_res
    echo $PWD
    sbatch $job_template_filename
    cd ..
  done 
  cd $top_level_dir
done
