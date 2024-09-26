# Set sources
# set(SOURCES 
#   strangsplitting.h
#   strangsplitting.cc)

# Define executable
add_executable(projects.cylinder_test.read_cylinder_mesh read_cylinder_mesh.cc ${SOURCES})

# link executable with required modules

target_link_libraries(projects.cylinder_test.read_cylinder_mesh
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

set_target_properties(projects.cylinder_test.read_cylinder_mesh        PROPERTIES FOLDER "projects/cylinder_test")

add_custom_target(projects.cylinder_test.read_cylinder_mesh_run         COMMAND splitting)

set_target_properties(projects.cylinder_test.read_cylinder_mesh_run         PROPERTIES FOLDER "projects/cylinder_test")

