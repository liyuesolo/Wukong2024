cmake_minimum_required(VERSION 3.20)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

project(WuKong2024 CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
        eigen #
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip #
        URL_HASH MD5=994092410ba29875184f7725e0371596 #
)
FetchContent_Populate(eigen)
add_library(eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS eigen)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/_deps/eigen-src/cmake/")
include_directories(${eigen_SOURCE_DIR})
add_compile_definitions(EIGEN_DONT_PARALLELIZE EIGEN_CHOLMOD_SUPPORT)

find_package (OpenMP REQUIRED COMPONENTS C CXX)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_subdirectory("Deps/polyscope")

include(libigl)
# igl_include(imgui)
# igl_include(glfw)
# igl_include(opengl)
igl_include(predicates)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(CHOLMOD REQUIRED)
if(CHOLMOD_FOUND)
   set(SPARSE_LIBS ${SPARSE_LIBS} ${CHOLMOD_LIBRARIES} ${EIGEN_BLAS_LIBRARIES} ${EIGEN_LAPACK_LIBRARIES})
   set(CHOLMOD_ALL_LIBS  ${CHOLMOD_LIBRARIES} ${EIGEN_BLAS_LIBRARIES} ${EIGEN_LAPACK_LIBRARIES})
endif()
include_directories(${CHOLMOD_INCLUDES})

add_definitions(-DEIGEN_USE_MKL_ALL)
add_definitions(-DMKL_LP64)


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${MKLROOT}/include)


add_subdirectory(Libs)
add_subdirectory(Projects)


message("**************************************************")
message("C++ Cmake Flags: ${CMAKE_CXX_FLAGS}")
message("**************************************************")
