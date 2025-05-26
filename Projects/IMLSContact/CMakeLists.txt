if (APPLE)
    set(DEPS boost_filesystem tbb 
        tbbmalloc tbbmalloc_proxy 
        igl::core
        suitesparseconfig metis cholmod amd camd ccolamd colamd
        gmp mpfr
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


add_executable(IMLSContact ${HEADERS} ${SOURCES})
add_library(IMLSContactLib STATIC ${HEADERS} ${SOURCES})


if (APPLE)
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    find_package(TBB REQUIRED)
    target_link_libraries(IMLSContact ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()

# find_package(Ceres REQUIRED)

target_include_directories(IMLSContact PUBLIC ../../Libs/spectra/include)
target_include_directories(IMLSContact PUBLIC ../../Deps/libigl/include)
target_include_directories(IMLSContact PUBLIC ../../Deps/geometry-central/include)
target_link_libraries(IMLSContact ${DEPS} polyscope geometry-central)