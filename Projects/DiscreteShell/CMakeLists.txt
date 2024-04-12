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




file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")


add_executable(DiscreteShell ${HEADERS} ${SOURCES})
add_library(DiscreteShellLib STATIC ${HEADERS} ${SOURCES})

if (APPLE)
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    find_package(TBB REQUIRED)
    target_link_libraries(DiscreteShell ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()

target_include_directories(DiscreteShell PUBLIC ../../Libs/spectra/include)
target_include_directories(DiscreteShell PUBLIC ../../Deps/libigl/include)
target_link_libraries(DiscreteShell ${DEPS} polyscope)


target_include_directories(DiscreteShellLib PUBLIC ../../Libs/spectra/include)
target_include_directories(DiscreteShellLib PUBLIC ../../Deps/libigl/include)
target_link_libraries(DiscreteShellLib ${DEPS} polyscope)
