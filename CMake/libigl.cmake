include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.4.0
)
FetchContent_MakeAvailable(libigl)
set(FETCHCONTENT_SOURCE_DIR_LIBIGL "${CMAKE_CURRENT_BINARY_DIR}/_deps/libigl")