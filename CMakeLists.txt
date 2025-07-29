cmake_minimum_required(VERSION 3.11)
project(SDTFE LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# PIC for everything (could be made target specific)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Sanitizer
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address -g")
# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer -fsanitize=address -g")
# set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address -g")

# Load QHULL
include(FetchContent)
FetchContent_Declare(
    qhull
    GIT_REPOSITORY https://github.com/qhull/qhull.git
    GIT_TAG v8.0.2
)
FetchContent_GetProperties(qhull)
if(NOT qhull_POPULATED)
    FetchContent_MakeAvailable(qhull)
    # patch incorrect cmake_source_dir in qhull
    execute_process(COMMAND sed -i "s@CMAKE_SOURCE_DIR@CMAKE_CURRENT_SOURCE_DIR@g" ${qhull_SOURCE_DIR}/CMakeLists.txt)
endif()

# load tiff
include(FetchContent)
FetchContent_Declare(
    libtiff
    GIT_REPOSITORY https://gitlab.com/libtiff/libtiff.git
    GIT_TAG v4.7.0
)
FetchContent_GetProperties(libtiff)
if(NOT libtiff_POPULATED)
    FetchContent_MakeAvailable(libtiff)
endif()

# OpenMP  (if found: OpenMP_FOUND)
find_package(OpenMP)

# set sources
set(sources
    src/bitmap.h src/bitmap.c
    src/constants.h
    src/dtfe.h src/dtfe.c
    src/io.h src/io.c
    src/triangulate.h src/triangulate.c
    src/utils.c src/utils.h
)

# library
add_library(libdtfe STATIC ${sources})
add_library(dtfe::libdtfe ALIAS libdtfe)
target_include_directories(libdtfe INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_include_directories(libdtfe PRIVATE ${qhull_SOURCE_DIR}/src/libqhull ${qhull_SOURCE_DIR}/src/libqhull_r)
target_link_libraries(libdtfe PRIVATE qhullstatic tiff ${M_LIBRARIES})
if(OpenMP_FOUND)
    target_link_libraries(libdtfe PRIVATE OpenMP::OpenMP_C)
endif()
set_target_properties(libdtfe PROPERTIES OUTPUT_NAME dtfe)

# executable
add_executable(dtfe src/main.c)
target_link_libraries(dtfe libdtfe)

# python library
add_subdirectory(python)
