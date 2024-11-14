add_library(alglib STATIC
    ap.cpp
    interpolation.cpp
    alglibinternal.cpp
    alglibmisc.cpp
    linalg.cpp
    solvers.cpp
    diffequations.cpp
    dataanalysis.cpp
    fasttransforms.cpp
    integration.cpp
    kernels_avx2.cpp
    kernels_fma.cpp
    kernels_sse2.cpp
    optimization.cpp
    specialfunctions.cpp
    statistics.cpp
)

target_include_directories(alglib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)