set(FETCHCONTENT_SOURCE_DIR_LIBIGL "${CMAKE_CURRENT_SOURCE_DIR}/../libigl")
include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)