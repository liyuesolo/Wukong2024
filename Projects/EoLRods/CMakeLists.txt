set(DEPS boost_filesystem tbb 
    tbbmalloc tbbmalloc_proxy 
    mkl_intel_lp64 
    mkl_sequential mkl_core 
    suitesparseconfig metis cholmod amd camd ccolamd colamd spqr
    gmp mpfr
    )

file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")

add_executable(EoLRods ${HEADERS} ${SOURCES})

target_include_directories(EoLRods PUBLIC ../../Libs/spectra/include)
target_link_libraries(EoLRods ${DEPS} polyscope)
