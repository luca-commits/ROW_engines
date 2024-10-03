# Set sources
# set(SOURCES 
#   strangsplitting.h
#   strangsplitting.cc)

set(SOURCES
 # read_cylinder_mesh.cc
  eddycurrent.h
  eddycurrent.cc 
  main_solve.cc
  utils.h
  )

# Define executable
add_executable(projects.cylinder_test.solve_eddy ${SOURCES})

# link executable with required modules

target_link_libraries(projects.cylinder_test.solve_eddy
  PUBLIC Boost::program_options  Eigen3::Eigen
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
  )

set_target_properties(projects.cylinder_test.solve_eddy         PROPERTIES FOLDER "projects/cylinder_test")

add_custom_target(projects.cylinder_test.solve_eddy_run         COMMAND splitting)

set_target_properties(projects.cylinder_test.solve_eddy_run     PROPERTIES FOLDER "projects/cylinder_test")

