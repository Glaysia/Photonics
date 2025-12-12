cmake_minimum_required(VERSION 3.28)
project(Photonics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_RPATH "$ENV{HOME}/.local/lib")
set(CMAKE_INSTALL_RPATH "$ENV{HOME}/.local/lib")

find_package(PkgConfig QUIET)
option(PHOTONICS_ENABLE_MEEP "Enable MEEP (FDTD) integration via pkg-config" ON)
set(PHOTONICS_MEEP_PKGCONFIG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/meep/install/lib/pkgconfig" CACHE PATH "Path to directory containing meep.pc (set after building/installing MEEP)")

add_library(photonics_lattice STATIC lattice_geometry.cpp)
target_include_directories(photonics_lattice PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(photonics_lattice PUBLIC ctlgeom)

add_library(photonics_monitor_suite STATIC monitor_suite.cpp)
target_include_directories(photonics_monitor_suite PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_q_analyzer STATIC q_analyzer.cpp)
target_include_directories(photonics_q_analyzer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_simulation_config STATIC simulation_config.cpp)
target_include_directories(photonics_simulation_config PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(Photonics main.cpp mylib.cpp)
target_link_libraries(Photonics PRIVATE
    photonics_lattice
    photonics_monitor_suite
    photonics_q_analyzer
    photonics_simulation_config
)

add_executable(lattice_geometry_test tests/lattice_geometry_test.cpp)
target_link_libraries(lattice_geometry_test PRIVATE photonics_lattice)

add_executable(monitor_suite_tests tests/monitor_suite_tests.cpp)
target_link_libraries(monitor_suite_tests PRIVATE photonics_monitor_suite)

add_executable(q_analyzer_test q_analyzer_test.cpp)
target_link_libraries(q_analyzer_test PRIVATE photonics_q_analyzer)

if (PHOTONICS_ENABLE_MEEP)
    if (PKG_CONFIG_FOUND)
        if (EXISTS "${PHOTONICS_MEEP_PKGCONFIG_DIR}/meep.pc")
            set(ENV{PKG_CONFIG_PATH} "${PHOTONICS_MEEP_PKGCONFIG_DIR}:$ENV{PKG_CONFIG_PATH}")
        endif()
        pkg_check_modules(MEEP QUIET meep)
        if (MEEP_FOUND)
            foreach(target
                photonics_lattice
                photonics_monitor_suite
                photonics_simulation_config
                Photonics
                lattice_geometry_test
                monitor_suite_tests
            )
                target_include_directories(${target} PRIVATE ${MEEP_INCLUDE_DIRS})
                target_link_libraries(${target} PRIVATE ${MEEP_LIBRARIES})
                target_compile_options(${target} PRIVATE ${MEEP_CFLAGS_OTHER})
                if (MEEP_LDFLAGS)
                    target_link_options(${target} PRIVATE ${MEEP_LDFLAGS})
                endif()
                if (MEEP_LDFLAGS_OTHER)
                    target_link_options(${target} PRIVATE ${MEEP_LDFLAGS_OTHER})
                endif()
            endforeach()
        else()
            message(STATUS "MEEP not found via pkg-config; build/install it under external/meep (e.g., prefix to external/meep/install) or set PHOTONICS_MEEP_PKGCONFIG_DIR.")
        endif()
    else()
        message(STATUS "pkg-config not found; cannot locate MEEP.")
    endif()
endif()
