# Set sources
# set(SOURCES 
#   strangsplitting.h
#   strangsplitting.cc)

# Define executable
add_executable(projects.team30a.read_t30_mesh read_t30_mesh.cc ${SOURCES})

# link executable with required modules

target_link_libraries(projects.team30a.read_t30_mesh
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

set_target_properties(projects.team30a.read_t30_mesh         PROPERTIES FOLDER "projects/team30a")

add_custom_target(projects.team30a.read_t30_mesh_run         COMMAND splitting)

set_target_properties(projects.team30a.read_t30_mesh_run         PROPERTIES FOLDER "projects/team30a")

