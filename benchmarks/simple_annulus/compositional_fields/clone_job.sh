#!/bin/bash

interpolation=('bilinear' 'cell_average')
grid_resolution =('2' '3' '4' '5' '6')
param_temp_filename='simple_annulus_particles.prm'
job_temp_file='./job_template.sh'
job_template_filename='job_template.sh'
aspect='/home/hlevinso/local/code/aspect_new/aspect/build/aspect'

top_level_dir=`pwd`

for interpolation_scheme in ${interpolation[@]}; do
  cd $interpolation_scheme
  for grid_res in ${grid_resolution[@]}; do
    cd $grid_res
    echo $PWD
    cp $job_temp_file .
    
    sed -i "s/SBATCH -J.*/SBATCH -J solkz_$((grid_res))x$((grid_res))/" $job_template_filename
    sed -i "s/SBATCH -e.*/SBATCH -e solkz_$((grid_res))x$((grid_res)).err/" $job_template_filename 
    sed -i "s/SBATCH -o.*/SBATCH -o solkz_$((grid_res))x$((grid_res)).out/" $job_template_filename 
    if [ $grid_res == '6' ]; then
      sed -i "s/SBATCH -N.*/SBATCH -N1/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n8/" $job_template_filename 
    elif [ $grid_res == '7' ]; then
      sed -i "s/SBATCH -N.*/SBATCH -N1/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n8/" $job_template_filename 
    elif [ $grid_res == '8' ]; then
      sed -i "s/SBATCH -N.*/SBATCH -N2/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n32/" $job_template_filename 
    elif [ $grid_res == '9' ]; then
      sed -i "s/SBATCH -N.*/SBATCH -N4/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n64/" $job_template_filename 
    elif [ $grid_res == '10' ]; then
      sed -i "s/SBATCH -N.*/SBATCH -N16/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n256/" $job_template_filename 
    else
      sed -i "s/SBATCH -N.*/SBATCH -N1/" $job_template_filename 
      sed -i "s/SBATCH -n.*/SBATCH -n1/" $job_template_filename 
    fi
    
    abs_path_of_param_file=`readlink -e $param_temp_filename` 
    echo "mpirun $aspect $abs_path_of_param_file" >> $job_template_filename

    cd ..
  done 
  cd $top_level_dir
done
