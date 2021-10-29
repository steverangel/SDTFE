cmake_minimum_required(VERSION 3.11)
project(SDTFE LANGUAGES C)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)


# Load QHULL
include(FetchContent)
FetchContent_Declare(
    qhull
    GIT_REPOSITORY https://github.com/qhull/qhull.git
    GIT_TAG v8.0.2
)
FetchContent_GetProperties(qhull)
if(NOT qhull_POPULATED)
    FetchContent_Populate(qhull)
    add_subdirectory(${qhull_SOURCE_DIR} ${qhull_BINARY_DIR})
endif()

# load tiff
include(FetchContent)
FetchContent_Declare(
    libtiff
    GIT_REPOSITORY https://gitlab.com/libtiff/libtiff.git
    GIT_TAG v4.3.0
)
FetchContent_GetProperties(libtiff)
if(NOT libtiff_POPULATED)
    FetchContent_Populate(libtiff)
    add_subdirectory(${libtiff_SOURCE_DIR} ${libtiff_BINARY_DIR})
endif()

# math
find_package(m)

# set sources
set(sources
    src/bitmap.h src/bitmap.c
    src/constants.h
    src/dtfe.h src/dtfe.c
    src/io.h src/io.c
    src/triangulate.h src/triangulate.c
    src/utils.c src/utils.h
    src/main.c
)

add_executable(dtfe ${sources})
target_include_directories(dtfe PRIVATE ${qhull_SOURCE_DIR}/src/libqhull)
target_link_libraries(dtfe PRIVATE qhullstatic tiff ${M_LIBRARIES})
