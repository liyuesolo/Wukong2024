set(DEPS boost_filesystem tbb 
    tbbmalloc tbbmalloc_proxy 
    mkl_intel_lp64 
    mkl_sequential mkl_core 
    suitesparseconfig metis cholmod amd camd ccolamd colamd spqr
    gmp mpfr
    )


include(FetchContent)
FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG v1.1.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc_toolkit)

file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")


add_executable(FEM3D ${HEADERS} ${SOURCES})

target_include_directories(FEM3D PUBLIC ../../Libs/spectra/include)
target_link_libraries(FEM3D ${DEPS} polyscope ipc::toolkit)
