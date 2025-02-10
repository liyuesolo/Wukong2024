if (APPLE)
    set(DEPS  
        igl::core
        )
else()
    set(DEPS boost_filesystem tbb 
        tbbmalloc tbbmalloc_proxy 
        mkl_intel_lp64 
        mkl_sequential mkl_core 
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
    )
endif()


file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")


add_executable(Homogenization ${HEADERS} ${SOURCES})

target_include_directories(Homogenization PUBLIC ../../Libs/spectra/include)
target_include_directories(Homogenization PUBLIC ../../Deps/libigl/include)
target_include_directories(Homogenization PUBLIC ../DiscreteShell/include)


if (APPLE)
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    find_package(TBB REQUIRED)
    target_link_libraries(Homogenization ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()

target_link_libraries(Homogenization ${DEPS} polyscope DiscreteShellLib)
