
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
ROW Wanner methods promise to avoid expensive iterations in timestepping, while still having a high order. This repository is used to put ROW on the test for the simulation of Eddy currents in electrical machines. ROW methods have been already used to simulate Eddy currents, but this has been limited to static meshes. Here we adapt the method to rotating geometries. 


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

## Usage

The configuration of a simulation is managed with inputs file. 
* The input of BDF-1/2 is: [input-iterative](projects/cylinder_test/xinput_newton.txt)
* The input of ROW static mesh is: [input-ROW-static](projects/cylinder_test/xinput.txt)
* The input of ROW rotating mesh is: [input-ROW-rotating](projects/cylinder_test/xinput_rotating.txt)
  

```sh
$ standard-readme-spec
# Prints out the standard-readme spec
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

