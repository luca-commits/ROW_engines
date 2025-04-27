
# Master Thesis: Adaptive timestepping for the Simulation of Electric Machines

This repository has been forked from the LehrFEM++ repository, all the files relevant to the Master thesis can be found in the folder [cylinder_test](projects/cylinder_test)

This repository contains:

1. The LehrFEM++ library
2. A program used to solve the Eddy current equation with iterative methods [BDF-1/2](projects/cylinder_test/solve_non-linear.cc)
3. Programs used to solve the Eddy current equation with ROW methods on a static mesh [ROW-static](projects/cylinder_test/solve_ROW_no_rotation_main.cc) and a rotating mesh [ROW-rotating](projects/cylinder_test/solve_ROW_complete.cc)
4. A script to compute the L2 norm of a solution computed with iterative methods or ROW [benchmarking](projects/cylinder_test/compute_L2_norm.sh). 

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)


## Background
ROW Wanner methods promise to avoid expensive iterations in timestepping, while still having a high order. This repository is used to put ROW to the test for the simulation of Eddy currents in electrical machines. ROW methods have been already used to simulate Eddy currents, but this has been limited to static meshes. Here we adapt the method to rotating geometries. 


## Install
First, LehrFEM++ needs to be installed, which can be done with CMake. After you clone the repository, create a build folder at the top level, and run CMake: 

```sh
$ clone https://github.com/luca-commits/ROW_engines
$ mkdir build
$ cd build
$ cmake ..
```
Then move the folder of this project and compile the code: 

```sh
$ cd projects/cylinder_test
$ make -j16
```

Then the folder structure needs to be reproduced in the build folder:

```
vtk_files/
├── time_dependent/
│   └── non-linear/
│       └── box_irregular_ring/
│           ├── bdf1/
│           └── bdf2/
│       └── transformator/
│           ├── bdf1/
│           ├── bdf2/
│       └── row_dynamic/
│       └── row_static/
```
which can be done with the following command: 
```
$ zsh create_folder_structure.sh
```

## Usage

The configuration of a simulation is managed with inputs file. 
* The input of BDF-1/2 is: [input-iterative](projects/cylinder_test/xinput_newton.txt)
* The input of ROW static mesh is: [input-ROW-static](projects/cylinder_test/xinput.txt)
* The input of ROW rotating mesh is: [input-ROW-rotating](projects/cylinder_test/xinput_rotating.txt)

For BDF-1/2 the input file looks like:
```
# step size
1e-3
# total time
2e-2       
# mesh, options: transformator, box_irregular_ring
transformator
# max_current
1e7
# conductivity     
5.96e7     
# exitation type, options: sinusoidal, ramp_up
sinusoidal
# exitation current parameter
50
# timestepping method    
bdf_2        
#geometry type, options: transformer, rotating_ring
transformer
#test mode
1
#adaptive // this parameter is only considere for Rosenbrock Wanner
1
#how frequently the timesteps are saved to vtk files
10. 
```
Here comes the part that is a bit confusing: 
* When using in test mode (creating the baseline, set test_mode = 1), step_size indicates the intervals at which the benchmark file is saved, and the actual step_size is computed by "how frequently the timesteps are saved to vtk files", so for example step_size=1e-3 and "how frequently the timesteps are saved to vtk files" = 10., results in step_size = 1e-4
* When using test_mode = 0 (actually running BDF-2 to see how it compares to ROW) step_size actually indicates the step size, and "how frequently the timesteps are saved to vtk files" should be set to match the timesteps at which the baseline has been recorded. For example if I have set the step_size to 1e-3 in baseline mode (test_mode=1), and I want to run BDF-2 with step_size 1e-4, then I need to set "how frequently the timesteps are saved to vtk files" to 10.
 
For the ROW input file, it is simpler, since the "saving step size" is explicitly defined: 

```
# step size
6e-5
# total time
2e-2       
# mesh, options: transformator, box_irregular_ring
transformator
# max_current
1e7
# conductivity     
5.96e7     
# exitation type, options: sinusoidal, ramp_up
sinusoidal
# exitation current parameter
50
# timestepping method, either bdf_1 or bdf_2
ROW        
transformer
#test mode
0
#adaptive // this parameter is only considere for Rosenbrock Wanner
0
#saving step size, how often benchmark and visualization files are saved
1e-3

```
For ROW, some arguments are ignored, e.g. timestepping method, and test_mode. 
Explenation of the arguments: 
* max_current: amplitude of the current, when using sinusoidal, an current to which the exitation linearly grows, when using ramp_up
* conductivity: in the geometry there are some conductive parts, like the conductive ring in the rotating geometry and the secondary coils in the transformer. This parameter sets the value of the conductivity in those regions.
* exitation type: ramp_up starts at 0, and linearly grows for a period defined by the exitation current parameter, and then stays constant.
* exitation parameter: defines frequency in sinusoidal, and growing time for ramp up
* adaptive: if 1 ROW will use adaptive timestepping, starting at step size defined by step_size parameter, and fixed timestep otherwise

## Running the programs

To run ROW on a static mesh: 
```
$ ./projects.cylinder_test.row_static
```
To run ROW on a rotating mesh:
```
$ ./projects.cylinder_test.row_dynamic
```
To run BDF-1/2:
```
$ ./projects.cylinder_test.solve_non-linear
```

[![Build Status](https://github.com/craffael/lehrfempp/workflows/Continuous%20Integration/badge.svg?branch=master)](https://github.com/craffael/lehrfempp/actions)

# LehrFEM++
Simple C++ Finite Element Framework for research and eduction optimzed for clarity and
flexibility with some trade-off concerning performance. This libary is used for the course _Numerical Methods for Partial Differential Euqations_ taught by Prof. R. Hiptmair at ETH Zurich.

* LehrFEM++ follows the [Google C++ Style
Guide](https://google.github.io/styleguide/cppguide.html#Naming).
* Adhere to the LehrFEM [coding style
  guidelines](https://github.com/craffael/lehrfempp/wiki/Contribute).
* Whenever adding core functionality, thorough testing is mandatory, following the
  instructions in the [testing
  guidelines](https://github.com/craffael/lehrfempp/wiki/Contribute).
* [Doxygen Class Documentation](https://craffael.github.io/lehrfempp)
* [Getting Started Guide](https://craffael.github.io/lehrfempp/getting_started.html)

## Contributors
- Raffael Casagrande (core developer)
- Ralf Hiptmair (core developer)
- Tobias Rohner (`projects/ipdg_stokes`, hp-fem in `lf::fe`)
- Anian Ruoss (Second order Geometry, Mesh Generators)
- Philippe Peter (`projects/dpg`)
- Amélie Justine Loher (`projects/FisherKPP`)
- Gina Magnussen (TIKZ output)
- Julien Gacon (`lf::base::comm`)
- Simon Meierhans

