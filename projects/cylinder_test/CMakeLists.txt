cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)
# Find BLAS
find_package(BLAS REQUIRED)

add_subdirectory(extern/alglib)

# Set sources
set(SOURCES
  # read_cylinder_mesh.cc
  eddycurrent.h
  eddycurrent.cc 
  implicit_euler.h
  utils.h
  BICG_stab.hpp
  remesh_airgap.hpp
  rotate_mesh.hpp
)

# Define executable for solving eddy current problem (steady-state)
# add_executable(projects.cylinder_test.solve_eddy main_solve.cc ${SOURCES})

# Define executable for solving eddy current problem (time-dependent)
# add_executable(projects.cylinder_test.solve_eddy_td timestepping_main.cc ${SOURCES})
# add_executable(projects.cylinder_test.read_mesh read_cylinder_mesh.cc)
# add_executable(projects.cylinder_test.solve_with_rotation solve_with_rotation.cc ${SOURCES})
add_executable(projects.cylinder_test.solve_non-linear solve_non-linear.cc ${SOURCES})

add_executable(projects.cylinder_test.test_meshes test_meshes.cpp  )
# add_executable(projects.cylinder_test.test_dof_mapping test_dof_mapping.cpp ${SOURCES})
add_executable(projects.cylinder_test.test_permeability test_permeability.cpp )

target_include_directories(projects.cylinder_test.test_meshes PUBLIC /usr/local/include/eigen3)
target_include_directories(projects.cylinder_test.test_meshes PUBLIC src/main)
target_include_directories(projects.cylinder_test.test_meshes PUBLIC /usr/local/include)

target_link_libraries(projects.cylinder_test.test_meshes /usr/local/lib/libgmsh.dylib)
target_link_libraries(projects.cylinder_test.test_meshes /usr/local/lib/libgmsh.dylib)

# target_include_directories(projects.cylinder_test.test_dof_mapping PUBLIC /usr/local/include)

target_link_libraries(projects.cylinder_test.test_permeability PRIVATE alglib)

# target_link_libraries(projects.cylinder_test.test_dof_mapping
#   PUBLIC Boost::program_options Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.uscalfe
#   BLAS::BLAS
#   /usr/local/lib/libgmsh.dylib
# )


# Link executables with required modules and BLAS
# target_link_libraries(projects.cylinder_test.read_mesh
#   PUBLIC Boost::program_options Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.uscalfe
#   BLAS::BLAS
# )

# Link executables with required modules and BLAS
# target_link_libraries(projects.cylinder_test.solve_eddy
#   PUBLIC Boost::program_options Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.uscalfe
#   BLAS::BLAS
# )

# target_link_libraries(projects.cylinder_test.solve_eddy_td
#   PUBLIC Boost::program_options Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.uscalfe
#   BLAS::BLAS
# )

# target_link_libraries(projects.cylinder_test.solve_with_rotation
#   PUBLIC Boost::program_options Eigen3::Eigen
#   lf.assemble
#   lf.base
#   lf.geometry
#   lf.io
#   lf.mesh
#   lf.mesh.hybrid2d
#   lf.mesh.test_utils
#   lf.mesh.utils
#   lf.refinement
#   lf.uscalfe
#   BLAS::BLAS
#   /usr/local/lib/libgmsh.dylib
# )

target_link_libraries(projects.cylinder_test.solve_non-linear
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
  /usr/local/lib/libgmsh.dylib
)

# Set target properties
# set_target_properties(projects.cylinder_test.solve_eddy            PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.solve_eddy_td         PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.solve_with_rotation   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.solve_non-linear   PROPERTIES FOLDER "projects/cylinder_test")

set_target_properties(projects.cylinder_test.test_meshes   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.test_dof_mapping   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.test_permeability   PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.test_dof_mapping   PROPERTIES FOLDER "projects/cylinder_test")

# target_include_directories(projects.cylinder_test.solve_with_rotation PUBLIC /usr/local/include)
target_include_directories(projects.cylinder_test.solve_non-linear PUBLIC /usr/local/include)

# Add custom targets for running
# add_custom_target(projects.cylinder_test.solve_eddy_run     COMMAND static)
# add_custom_target(projects.cylinder_test.solve_eddy_td_run  COMMAND static)
# add_custom_target(projects.cylinder_test.solve_with_rotation_run  COMMAND static)
# add_custom_target(projects.cylinder_test.test_meshes_run  COMMAND static)
# add_custom_target(projects.cylinder_test.test_dof_mapping_run  COMMAND static)
add_custom_target(projects.cylinder_test.solve_non_linear_run  COMMAND static)

# Set folder properties for the custom targets
# set_target_properties(projects.cylinder_test.solve_eddy_run        PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.solve_eddy_td_run     PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.solve_with_rotation_run     PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.test_meshes     PROPERTIES FOLDER "projects/cylinder_test")
# set_target_properties(projects.cylinder_test.test_dof_mapping_run     PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.solve_non_linear_run    PROPERTIES FOLDER "projects/cylinder_test")

