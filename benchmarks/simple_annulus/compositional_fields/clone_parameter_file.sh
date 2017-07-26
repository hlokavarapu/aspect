#!/bin/bash

#The various interpolation schemes and grid resolution tests.
# Note that, one could add or delete to grid_resolution to run various models. 
interpolation=('bilinear' 'cell_average')
grid_resolution=('2', '3', '4', '5', '6')
#list_of_number_of_particles_per_cell_in_x=('16')
#list_of_number_of_particles_per_cell_in_y=('2','4','8','16','32','64','128')
#list_of_number_of_particles_per_cell_in_total=($((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)),$((list_of_number_of_particles_per_cell_in_x*list_of_number_of_particles_per_cell_in_y)))

#These variables are user specific and will need to be appropriately changed:
# 
# Assuming that setup of models will be done from the same working directory as the script is executed.
top_level_dir=`echo $PWD`
#The name of the template parameter file.
param_temp_filename="simple_annulus_particles.prm"
#Assuming that the parameter file is in the present working directory...
parameter_template_file=`echo ${PWD}/${param_temp_filename}`
# Change this to the absolute path to compiled benchmark library. 
# Make sure that each forward slash is delimted by a back slash. I.e. '\/'
path_to_benchmark_library='./libsimple_annulus_compositional_fields.so'

for interpolation_scheme in ${interpolation[@]}; do
  if [ ! -d $interpolation_scheme ]; then
    mkdir $interpolation_scheme
  fi

  cd $interpolation_scheme
  for grid_res in ${grid_resolution[@]}; do
    if [ ! -d $grid_res ]; then
      mkdir $grid_res 
    fi 
    second_level_dir=`echo $PWD`
    cd $grid_res
    for n_p_x in ${list_of_number_of_particles_per_cell_in_x[@]}; do
      mkdir $n_p_x
      cd $n_p_x
      echo $PWD
      
      cp ${parameter_template_file} ./${param_temp_filename}

      sed -i "s/set Additional shared libraries.*/set Additional shared libraries = ${path_to_benchmark_library}/" $param_temp_filename
      sed -i "s/set Output directory.*/set Output directory = $((grid_res))x$((grid_res))_$((n_p_x))_output/" $param_temp_filename  
      sed -i "s/set Initial global refinement.*/set Initial global refinement = $((grid_res))/" $param_temp_filename      
       
       if [ $interpolation_scheme == 'cell_average' ]; then
         sed -i "s/set Interpolation scheme.*/set Interpolation scheme = cell average/" $param_temp_filename
       elif [ $interpolation_scheme == 'bilinear' ]; then
         sed -i "s/set Interpolation scheme.*/set Interpolation scheme = bilinear/" $param_temp_filename
       elif [ $interpolation_scheme == 'biquadratic' ]; then
         sed -i "s/set Interpolation scheme.*/set Interpolation scheme = biquadratic/" $param_temp_filename
       fi
   
       min_x=`echo "1/($grid_res*$n_p_x*2)" | bc -l`
       max_x=`echo "1 - 1/($grid_res*$n_p_x*2)" | bc -l`
       min_y=`echo "1/($grid_res*$n_p_x*2)" | bc -l`
       max_y=`echo "1 - 1/($grid_res*$n_p_x*2)" | bc -l`
     
       #sed -i "s/set Number of tracers.*/set Number of tracers = $((grid_res*grid_res*n_p_x*n_p_x))/" $param_temp_filename
       #sed -i "s/set Number of particles in x.*/set Number of particles in x = $((grid_res*n_p_x))/" $param_temp_filename
      # sed -i "s/set Number of particles in y.*/set Number of particles in y = $((grid_res*n_p_x))/" $param_temp_filename
     
       delta_x=`echo "1/($grid_res * $n_p_x)" | bc -l`
       delta_y=`echo "1/($grid_res * $n_p_x)" | bc -l`
     
       sed -i "s/set Delta x.*/set Delta x = $delta_x/" $param_temp_filename
       sed -i "s/set Delta y.*/set Delta y = $delta_y/" $param_temp_filename
       
       cd ..
    done
    cd $second_level_dir
  done 
  cd $top_level_dir
done
