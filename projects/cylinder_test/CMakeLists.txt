# Find BLAS
find_package(BLAS REQUIRED)

# Set sources
set(SOURCES
  # read_cylinder_mesh.cc
  eddycurrent.h
  eddycurrent.cc 
  eddycurrent_time_dependent.h
  eddycurrent_time_dependent.cc
  implicit_euler.h
  utils.h
)

# Define executable for solving eddy current problem (steady-state)
add_executable(projects.cylinder_test.solve_eddy main_solve.cc ${SOURCES})

# Define executable for solving eddy current problem (time-dependent)
add_executable(projects.cylinder_test.solve_eddy_td timestepping_main.cc ${SOURCES})

# Link executables with required modules and BLAS
target_link_libraries(projects.cylinder_test.solve_eddy
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  BLAS::BLAS
)

target_link_libraries(projects.cylinder_test.solve_eddy_td
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  BLAS::BLAS
)

# Set target properties
set_target_properties(projects.cylinder_test.solve_eddy            PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.solve_eddy_td         PROPERTIES FOLDER "projects/cylinder_test")

# Add custom targets for running
add_custom_target(projects.cylinder_test.solve_eddy_run     COMMAND static)
add_custom_target(projects.cylinder_test.solve_eddy_td_run  COMMAND static)

# Set folder properties for the custom targets
set_target_properties(projects.cylinder_test.solve_eddy_run        PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.solve_eddy_td_run     PROPERTIES FOLDER "projects/cylinder_test")
