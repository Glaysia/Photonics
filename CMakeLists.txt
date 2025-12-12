cmake_minimum_required(VERSION 3.28)
project(Photonics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_RPATH "$ENV{HOME}/.local/lib")
set(CMAKE_INSTALL_RPATH "$ENV{HOME}/.local/lib")

find_package(PkgConfig QUIET)
option(PHOTONICS_ENABLE_MEEP "Enable MEEP (FDTD) integration via pkg-config" ON)
set(PHOTONICS_MEEP_PKGCONFIG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/meep/install/lib/pkgconfig" CACHE PATH "Path to directory containing meep.pc (set after building/installing MEEP)")

add_executable(Photonics
    main.cpp
    mylib.cpp
    simulation_config.cpp
)

target_include_directories(Photonics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (PHOTONICS_ENABLE_MEEP)
    if (PKG_CONFIG_FOUND)
        if (EXISTS "${PHOTONICS_MEEP_PKGCONFIG_DIR}/meep.pc")
            set(ENV{PKG_CONFIG_PATH} "${PHOTONICS_MEEP_PKGCONFIG_DIR}:$ENV{PKG_CONFIG_PATH}")
        endif()
        pkg_check_modules(MEEP QUIET meep)
        if (MEEP_FOUND)
            target_include_directories(Photonics PRIVATE ${MEEP_INCLUDE_DIRS})
            target_link_libraries(Photonics PRIVATE ${MEEP_LIBRARIES})
            target_compile_options(Photonics PRIVATE ${MEEP_CFLAGS_OTHER})
            if (MEEP_LDFLAGS)
                target_link_options(Photonics PRIVATE ${MEEP_LDFLAGS})
            endif()
            if (MEEP_LDFLAGS_OTHER)
                target_link_options(Photonics PRIVATE ${MEEP_LDFLAGS_OTHER})
            endif()
        else()
            message(STATUS "MEEP not found via pkg-config; build/install it under external/meep (e.g., prefix to external/meep/install) or set PHOTONICS_MEEP_PKGCONFIG_DIR.")
        endif()
    else()
        message(STATUS "pkg-config not found; cannot locate MEEP.")
    endif()
endif()
