#[=====[

CMake script for building dependencies of `rest_libcint`.

- Specify `CINT_DIR` as cmake argument or `REST_CINT_DIR` (or `REST_EXT_DIR`) as environment variable to find the library `libcint.so`.
- Specify `CINT_SRC` as cmake argument or `REST_CINT_SRC` as environment variable to source code of libcint.

]=====]

cmake_minimum_required(VERSION 3.15)
project(rest_dep)

include(ExternalProject)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

# == function to parse variables ==

macro(parsevar ARGNAME ARG_LIB ARG_DIR ARG_SRC ENV_DIR ENV_SRC)
    # 1. Try to find library if cache variables (-DARG_DIR=xxx, -DARG_SRC=xxx) are set
    set(${ARG_DIR} "" CACHE PATH "User specified directory of library (${ARGNAME})")
    set(${ARG_SRC} "" CACHE PATH "User specified directory of source code (${ARGNAME})")
    set(${ARG_LIB} "")
    # 2. If above variables are not defined, try to find library if environment variables are set
    #    (cache variables will override environment variables anyway)
    if(NOT ${ARG_SRC} AND (NOT "$ENV{${ENV_SRC}}" EQUAL ""))
        set(${ARG_SRC} $ENV{${ENV_SRC}})
    endif()
    # 3. Find by library if directory is specified
    if(${ARG_DIR} OR (NOT "$ENV{${ENV_DIR}}" EQUAL ""))
        message(STATUS "User specified directory of library (${ARGNAME}): ${${ARG_DIR}}")
        message(STATUS "Environmental specified directory of library (${ARGNAME}): $ENV{${ENV_DIR}}")
        find_library(LIB_FOUND NAMES ${ARGNAME} HINTS "${${ARG_DIR}}" "$ENV{${ENV_DIR}}" "$ENV{REST_EXT_DIR}")
        if(LIB_FOUND)
            message(STATUS "Found library (${ARGNAME}): ${LIB_FOUND}")
            set(${ARG_LIB} ${LIB_FOUND})
        else()
            message(STATUS "Not found library (${ARGNAME}), try to build from source.")
        endif()
    endif()
    # After this function
    # - ${ARG_LIB} will be set if library is found
    # - ${ARG_SRC} will be set if source code is specified
endmacro()

# == libcint ==

parsevar("libcint.so" CINT_LIB CINT_DIR CINT_SRC REST_CINT_DIR REST_CINT_SRC)
set(CINT_VER "v6.1.2" CACHE STRING "Version of library (libcint.so)")

# build with git source
if(NOT CINT_LIB)
    if(CINT_SRC)
        message(STATUS "User specified directory of source code (libcint.so): ${CINT_SRC}")
    else()
        set(CINT_SRC https://github.com/sunqm/libcint.git)
        message(STATUS "Download source code (libcint.so): ${CINT_SRC}")
    endif()

    ExternalProject_Add(libcint
        GIT_REPOSITORY ${CINT_SRC}
        GIT_TAG v6.1.2
        PREFIX ${PROJECT_BINARY_DIR}/deps
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=RELEASE
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DWITH_RANGE_COULOMB:STRING=1
            -DWITH_FORTRAN:STRING=0
            -DWITH_CINT2_INTERFACE:STRING=1
            -DWITH_F12:STRING=1
            -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
            # -DENABLE_STATIC=1
    )
    
    set(CINT_LIB ${PROJECT_BINARY_DIR}/deps/lib/libcint.so)
endif()

install(FILES ${CINT_LIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

