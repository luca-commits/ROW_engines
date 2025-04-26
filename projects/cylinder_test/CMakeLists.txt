cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)
# Find BLAS
find_package(BLAS REQUIRED)

add_subdirectory(extern/alglib)

# Set sources
set(SOURCES
  eddycurrent.h
  eddycurrent.cc 
  implicit_euler.h
  utils.h
  BICG_stab.hpp
  remesh_airgap.hpp
  rotate_mesh.hpp
  rosenbrock_wanner.hpp
  increments_no_rotations.hpp
)

add_executable(projects.cylinder_test.solve_non-linear solve_non-linear.cc ${SOURCES})
add_executable(projects.cylinder_test.row_static solve_ROW_no_rotation_main.cc ${SOURCES})
add_executable(projects.cylinder_test.row_dynamic solve_ROW_complete.cc ${SOURCES})
target_include_directories(projects.cylinder_test.solve_non-linear PUBLIC /usr/local/include)
target_include_directories(projects.cylinder_test.row_static PUBLIC /usr/local/include)
target_include_directories(projects.cylinder_test.row_dynamic PUBLIC /usr/local/include)


set_target_properties(projects.cylinder_test.solve_non-linear   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.row_static   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.row_dynamic   PROPERTIES FOLDER "projects/cylinder_test")



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
  PRIVATE alglib
)

target_link_libraries(projects.cylinder_test.row_static
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
  PRIVATE alglib
)

target_link_libraries(projects.cylinder_test.row_dynamic
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
  PRIVATE alglib
)


set_target_properties(projects.cylinder_test.row_static   PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.row_dynamic   PROPERTIES FOLDER "projects/cylinder_test")
add_custom_target(projects.cylinder_test.row_static_run  COMMAND static)

set_target_properties(projects.cylinder_test.row_static     PROPERTIES FOLDER "projects/cylinder_test")
set_target_properties(projects.cylinder_test.row_dynamic     PROPERTIES FOLDER "projects/cylinder_test")
