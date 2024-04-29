cmake_minimum_required(VERSION 3.20)
add_definitions(-w)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
endif()

project(WuKong2024 CXX)

message(STATUS "${CMAKE_BUILD_TYPE} Build")


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -O0 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O0 -fopenmp -Wno-unused-variable -Wno-unused-but-set-variable -Wsign-compare -Wreorder")


find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (APPLE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O0 -Wno-unused-variable -Wno-unused-but-set-variable -Wsign-compare -Wreorder")
else()
    find_package (OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -O0 -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O0 -fopenmp -Wno-unused-variable -Wno-unused-but-set-variable -Wsign-compare -Wreorder")
endif()


add_subdirectory("Deps/polyscope")
add_subdirectory("Deps/libigl")

include(libigl)
igl_include(predicates)

find_package(CGAL REQUIRED COMPONENTS Core)
  
if (NOT APPLE)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    add_definitions(-DMKL_LP64)

    set(MKLROOT "/opt/intel/oneapi/mkl/2024.1")
    include_directories(${MKLROOT}/include)
    link_directories(${MKLROOT}/lib)
endif()

include_directories(/usr/local/include/suitesparse)
include_directories(/usr/local/include/)
link_directories(/usr/local/lib)

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(Libs)
add_subdirectory(Projects)


message("**************************************************")
message("C++ Cmake Flags: ${CMAKE_CXX_FLAGS}")
message("**************************************************")
