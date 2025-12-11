cmake_minimum_required(VERSION 3.28)
project(Photonics)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Photonics main.cpp mylib.cpp)

target_include_directories(Photonics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Photonics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/mfem/install/include)

set(MFEM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/mfem/install/lib/cmake/mfem" CACHE PATH "Path to the MFEM package configuration")
find_package(mfem CONFIG REQUIRED)
target_link_libraries(Photonics PRIVATE mfem)
