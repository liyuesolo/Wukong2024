if (APPLE)
    set(DEPS  
        igl::core
        )
else()
    set(DEPS boost_filesystem tbb 
        tbbmalloc tbbmalloc_proxy 
        mkl_intel_lp64 igl::core
        mkl_sequential mkl_core 
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
        )
endif()

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

igl_include(copyleft tetgen)


add_executable(IMLSContact ${HEADERS} ${SOURCES})

target_include_directories(IMLSContact PUBLIC ../../Libs/spectra/include)
target_link_libraries(IMLSContact ${DEPS} polyscope ipc::toolkit igl_copyleft::tetgen)

find_package(CHOLMOD REQUIRED)
include_directories(${CHOLMOD_INCLUDES})
target_link_libraries(IMLSContact ${CHOLMOD_LIBRARIES})
