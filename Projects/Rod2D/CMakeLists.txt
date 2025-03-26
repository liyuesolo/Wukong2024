if(APPLE)
	set(DEPS igl::core)
else()
	set(DEPS boost_filesystem tbb 
             tbbmalloc tbbmalloc_proxy 
             mkl_intel_lp64 igl::core 
             mkl_sequential mkl_core  
             suitesparseconfig metis cholmod amd camd ccolamd colamd spqr 
             gmp mpfr)

endif()

file(GLOB HEADERS "include/*.h" "autodiff/*.h")
file(GLOB SOURCES "src/*.cpp" "autodiff/*.cpp")

add_executable(Rod2D ${HEADERS} ${SOURCES})
if(APPLE)
     find_package(CHOLMOD REQUIRED) 
     include_directories(${CHOLMOD_INCLUDES}) 
     find_package(TBB REQUIRED) 
	target_link_libraries(Rod2D ${CHOLMOD_LIBRARIES} TBB::tbb)
endif()
target_include_directories(Rod2D PUBLIC ../../Libs/spectra/include)
target_include_directories(Rod2D PUBLIC ../../Deps/libigl/include)
target_link_libraries(Rod2D ${DEPS} polyscope)