#!/bin/bash

# Set the base directory (current directory is /projects/cylinder_test/)
base_dir=$(pwd)

# Define the parent projects directory (../../build/projects/cylinder_test)
projects_dir="$base_dir/../../build/projects/cylinder_test/vtk_files/time_dependent/"

# Define the directories to create (based on the structure from previous answer)
directories=(
    "$projects_dir/non-linear"
    "$projects_dir/non-linear/box_irregular_ring"
    "$projects_dir/non-linear/box_irregular_ring/bdf1"
    "$projects_dir/non-linear/box_irregular_ring/bdf2"
    "$projects_dir/non-linear/transformator/bdf1"
    "$projects_dir/non-linear/transformator/bdf2"
    "$projects_dir/row_dynamic"
    "$projects_dir/row_static"
)

# Create the directories if they don't exist
for dir in "${directories[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
    echo "Created directory: $dir"
  else
    echo "Directory already exists: $dir"
  fi
done

echo "Folder structure replicated successfully under $projects_dir."