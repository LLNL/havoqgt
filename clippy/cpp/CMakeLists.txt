# Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy
# Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: MIT

# Works with 3.11 and tested through 3.15 (not tested yet)
cmake_minimum_required(VERSION 3.11...3.15)

project(CLIPPy
    VERSION 0.1
    DESCRIPTION "Command Line Interface Plus Python"
    LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    
    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    include(CTest)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND)
	#add_subdirectory(docs)
    else()
	message(STATUS "Doxygen not found, not building docs")
    endif()
endif()


#
#  Nlohmann JSON
#
include(ExternalProject)
ExternalProject_Add(nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.9.1/json.hpp
    DOWNLOAD_NO_EXTRACT 1
    LOG_DOWNLOAD 1
    DOWNLOAD_DIR nlohmann_json/include/nlohmann/
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    CONFIGURE_COMMAND ""
    )
set(NLOHMANN_JSON_INCLUDE_DIR ${PROJECT_BINARY_DIR}/nlohmann_json/include)

#
#  Metall
#
ExternalProject_Add(metall
    GIT_REPOSITORY https://github.com/LLNL/metall.git
    GIT_TAG v0.5
    CONFIGURE_COMMAND ""
    GIT_PROGRESS 1
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    #LOG_DOWNLOAD 1
    )
ExternalProject_Get_Property(metall source_dir)
set(METALL_INCLUDE_DIR ${source_dir}/include)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")

option(TEST_WITH_SLURM "Run tests with Slurm" OFF)

# Header-only library, so likely not have src dir 
# add_subdirectory(src)

# Testing & examples are only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
    #add_subdirectory(test)
    # Example codes are here.
    add_subdirectory(examples)
    #add_subdirectory(maptrix)
endif()
