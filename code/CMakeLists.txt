cmake_minimum_required(VERSION 3.7)

# set the project name and version
project(AletasFT CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# compiler flags that applies to all
if(MSVC)
  add_compile_options(/W4 /O2)
else()
  add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Wno-error=maybe-uninitialized
    -Wno-error=unknown-pragmas
    -O3
    -funroll-loops
    -march=native
    -mfpmath=sse)
endif()

add_compile_definitions(NDEBUG)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ThirdParties)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/include)
