set(CMAKE_CXX_STANDARD 17)

# load pybind11
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

# add module
pybind11_add_module(pydtfe pydtfe.cpp)
target_link_libraries(pydtfe PRIVATE libdtfe)
