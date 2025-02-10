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

add_executable(MMCCells ${HEADERS} ${SOURCES})

if (APPLE)
    find_package(CHOLMOD REQUIRED)
    include_directories(${CHOLMOD_INCLUDES})
    find_package(TBB REQUIRED)
    target_link_libraries(MMCCells ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()

target_include_directories(MMCCells PUBLIC ../../Libs/spectra/include)
target_link_libraries(MMCCells ${DEPS} polyscope)
