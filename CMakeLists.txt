cmake_minimum_required(VERSION 3.28)
project(Photonics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Photonics main.cpp mylib.cpp)

target_include_directories(Photonics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
