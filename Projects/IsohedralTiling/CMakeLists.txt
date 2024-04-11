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
    
file(GLOB HEADERS "include/*.h" "autodiff/*.h" "tactile/.*h" "tactile/.*hpp" )
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp" "tactile/*.cpp")

add_executable(IsohedralTiling ${HEADERS} ${SOURCES})

# target_include_directories(IsohedralTiling PUBLIC ../../Libs/spectra/include)
# target_include_directories(IsohedralTiling PUBLIC ../../Deps/libigl/include)

# add clipper library
add_library(clipper STATIC clipper/clipper.cpp)
target_include_directories(clipper PUBLIC clipper/clipper.h)
target_link_libraries(clipper ${DEPS})
set(DEPS ${DEPS} clipper)

# target_link_libraries(IsohedralTiling ${DEPS} polyscope)
# find_package(TBB REQUIRED)
# target_link_libraries(IsohedralTiling TBB::tbb)

find_package(Gmsh REQUIRED)
include_directories(${GMSH_INC})
target_link_libraries(IsohedralTiling ${GMSH_LIB})
