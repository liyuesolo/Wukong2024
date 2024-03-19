set(DEPS boost_filesystem tbb 
    tbbmalloc tbbmalloc_proxy 
    mkl_intel_lp64 igl::core
    mkl_sequential mkl_core 
    suitesparseconfig metis cholmod amd camd ccolamd colamd spqr
    gmp mpfr
    )


file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")


add_executable(Homogenization ${HEADERS} ${SOURCES})

target_include_directories(Homogenization PUBLIC ../../Libs/spectra/include)
target_include_directories(Homogenization PUBLIC ../../Deps/libigl/include)
target_include_directories(Homogenization PUBLIC ../DiscreteShell/include)

target_link_libraries(Homogenization ${DEPS} polyscope DiscreteShellLib)
