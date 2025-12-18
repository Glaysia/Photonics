cmake_minimum_required(VERSION 3.28)
project(Photonics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Allow overriding the prefix where local dependencies (libctl, libharminv, etc.)
# have been installed.  The defaults match the build instructions under external/.
set(PHOTONICS_LOCAL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/external/local"
    CACHE PATH "Prefix containing locally built dependencies (libctl, harminv, ...)")
set(PHOTONICS_LOCAL_INCLUDE_DIR "${PHOTONICS_LOCAL_PREFIX}/include")
set(PHOTONICS_LOCAL_LIB_DIR "${PHOTONICS_LOCAL_PREFIX}/lib")
set(PHOTONICS_LOCAL_LIB64_DIR "${PHOTONICS_LOCAL_PREFIX}/lib64")

# Ensure libctl headers/libraries (and other local deps) are visible to every target.
if (PHOTONICS_LOCAL_INCLUDE_DIR)
    include_directories(${PHOTONICS_LOCAL_INCLUDE_DIR})
endif()
if (PHOTONICS_LOCAL_LIB_DIR OR PHOTONICS_LOCAL_LIB64_DIR)
    link_directories(${PHOTONICS_LOCAL_LIB_DIR} ${PHOTONICS_LOCAL_LIB64_DIR})
endif()

# Keep user-local libs available, but ensure project-provided MEEP (if used)
# takes precedence at runtime.
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

add_library(photonics_harminv_runner STATIC harminv_runner.cpp)
target_include_directories(photonics_harminv_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(photonics_harminv_runner PUBLIC photonics_q_analyzer)

add_library(photonics_simulation_config STATIC simulation_config.cpp)
target_include_directories(photonics_simulation_config PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_source_manager STATIC source_manager.cpp)
target_include_directories(photonics_source_manager PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_sweep_runner STATIC sweep_runner.cpp)
target_include_directories(photonics_sweep_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_geometry_export STATIC geometry_export.cpp)
target_include_directories(photonics_geometry_export PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(photonics_geometry_export PUBLIC photonics_lattice)

add_library(photonics_field_slice_export STATIC field_slice_export.cpp)
target_include_directories(photonics_field_slice_export PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(photonics_simulation_runner STATIC simulation_runner.cpp)
target_include_directories(photonics_simulation_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(photonics_simulation_runner PUBLIC photonics_simulation_config photonics_source_manager)

add_executable(Photonics main.cpp mylib.cpp)
target_link_libraries(Photonics PRIVATE
    photonics_lattice
    photonics_monitor_suite
    photonics_q_analyzer
    photonics_harminv_runner
    photonics_simulation_config
    photonics_source_manager
    photonics_sweep_runner
    photonics_geometry_export
    photonics_field_slice_export
    photonics_simulation_runner
)

add_executable(lattice_geometry_test tests/lattice_geometry_test.cpp)
target_link_libraries(lattice_geometry_test PRIVATE photonics_lattice)

add_executable(monitor_suite_tests tests/monitor_suite_tests.cpp)
target_link_libraries(monitor_suite_tests PRIVATE photonics_monitor_suite)

add_executable(q_analyzer_test q_analyzer_test.cpp)
target_link_libraries(q_analyzer_test PRIVATE photonics_q_analyzer)

add_executable(source_manager_tests source_manager_test.cpp)
target_link_libraries(source_manager_tests PRIVATE photonics_source_manager)

if (PHOTONICS_ENABLE_MEEP)
    if (PKG_CONFIG_FOUND)
        if (EXISTS "${PHOTONICS_MEEP_PKGCONFIG_DIR}/meep.pc")
            set(ENV{PKG_CONFIG_PATH} "${PHOTONICS_MEEP_PKGCONFIG_DIR}:$ENV{PKG_CONFIG_PATH}")
        endif()
        pkg_check_modules(MEEP QUIET meep)
        if (NOT MEEP_FOUND)
            set(PHOTONICS_MEEP_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/external/meep/install")
            if (EXISTS "${PHOTONICS_MEEP_PREFIX}/include/meep.hpp" AND
                (EXISTS "${PHOTONICS_MEEP_PREFIX}/lib/libmeep.so" OR EXISTS "${PHOTONICS_MEEP_PREFIX}/lib/libmeep.a"))
                set(MEEP_FOUND TRUE)
                set(MEEP_INCLUDE_DIRS "${PHOTONICS_MEEP_PREFIX}/include")
                set(MEEP_LIBRARIES "${PHOTONICS_MEEP_PREFIX}/lib/libmeep.so")
                set(MEEP_LIBRARY_DIRS "${PHOTONICS_MEEP_PREFIX}/lib")
                message(STATUS "Using bundled MEEP from ${PHOTONICS_MEEP_PREFIX} (pkg-config module not found).")
            endif()
        endif()
        if (MEEP_FOUND)
            foreach(target
                photonics_lattice
                photonics_monitor_suite
                photonics_harminv_runner
                photonics_simulation_config
                photonics_source_manager
                photonics_geometry_export
                photonics_field_slice_export
                photonics_simulation_runner
                Photonics
                lattice_geometry_test
                monitor_suite_tests
                source_manager_tests
            )
                target_include_directories(${target} PRIVATE ${MEEP_INCLUDE_DIRS})
                target_link_libraries(${target} PRIVATE ${MEEP_LIBRARIES})
                target_compile_options(${target} PRIVATE ${MEEP_CFLAGS_OTHER})
                if (MEEP_LIBRARY_DIRS)
                    # Prefer the selected MEEP libdir over ~/.local/lib so we don't
                    # accidentally pick up a different libmeep at runtime.
                    set_target_properties(${target} PROPERTIES
                        BUILD_RPATH "${MEEP_LIBRARY_DIRS};$ENV{HOME}/.local/lib"
                        INSTALL_RPATH "${MEEP_LIBRARY_DIRS};$ENV{HOME}/.local/lib"
                    )
                endif()
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
