cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)
project(rotate_mesh)

set(sources 
    rotate_mesh.hpp
    remesh_airgap.hpp
    )

###############################################################################
## target definitions #########################################################
###############################################################################

# add the data to the target, so it becomes visible in some IDE
add_executable(rotate ${sources} main.cpp)
add_executable(test_meshes ${sources} test_meshes.cpp)
add_executable(test_boundary_matching ${sources} test_boundary_matching.cpp)


# just for example add some compiler flags
target_compile_options(rotate PUBLIC -std=c++1y -Wall -Wfloat-conversion)
target_link_libraries(test_meshes /usr/local/lib/libgmsh.dylib)
target_link_libraries(test_boundary_matching /usr/local/lib/libgmsh.dylib)

# this lets me include files relative to the root source directory with a <> pair
target_include_directories(rotate PUBLIC src/main)
target_include_directories(test_meshes PUBLIC src/main)
target_include_directories(test_boundary_matching PUBLIC src/main)

# this copies all resource files in the build directory
# we need this, because we want to work with paths relative to the executable
file(COPY ${data} DESTINATION resources)

###############################################################################
## dependencies ###############################################################
###############################################################################

# this defines the variables Boost_LIBRARIES that contain all library names
# that we need to link to
# find_package(Boost 1.36.0 COMPONENTS filesystem system REQUIRED)
# set(Gmsh_DIR "/opt/homebrew/opt/gmsh")

# target_link_libraries(example PUBLIC
#   ${Boost_LIBRARIES}
#   # here you can add any library dependencies
# )

# Specify the Gmsh include and library paths manually
target_include_directories(rotate PUBLIC /usr/local/include)
target_include_directories(rotate PUBLIC /usr/local/include/eigen3)
target_include_directories(test_meshes PUBLIC /usr/local/include)
target_include_directories(test_meshes PUBLIC /usr/local/include/eigen3)
target_include_directories(test_boundary_matching PUBLIC /usr/local/include)
target_include_directories(test_boundary_matching PUBLIC /usr/local/include/eigen3)


target_link_libraries(rotate /usr/local/lib/libgmsh.dylib)
target_link_libraries(test_meshes /usr/local/lib/libgmsh.dylib)
target_link_libraries(test_boundary_matching /usr/local/lib/libgmsh.dylib)
