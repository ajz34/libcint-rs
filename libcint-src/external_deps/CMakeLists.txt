cmake_minimum_required(VERSION 3.15)
project(dep-libcint)

include(ExternalProject)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

set(CINT_SRC "" CACHE PATH "User specified directory of source code (libcint)")
if(NOT CINT_SRC AND (NOT $ENV{CINT_SRC} EQUAL ""))
    set(CINT_SRC $ENV{CINT_SRC})
endif()

if(CINT_SRC)
    message(STATUS "User specified directory of source code (libcint): ${CINT_SRC}")
else()
    # set(CINT_SRC https://github.com/sunqm/libcint.git)
    set(CINT_SRC /home/a/Git-Others/libcint)
    message(STATUS "Download source code (libcint): ${CINT_SRC}")
endif()

set(CINT_VER "" CACHE PATH "User specified version of library (git tag) (libcint)")
if(NOT CINT_VER AND (NOT $ENV{CINT_VER} EQUAL ""))
    set(CINT_VER $ENV{CINT_VER})
endif()

if(CINT_VER)
    message(STATUS "User specified version of library (libcint): ${CINT_VER}")
else()
    set(CINT_VER v6.1.2)
    message(STATUS "Download version of library (libcint): ${CINT_VER}")
endif()

ExternalProject_Add(cint
    GIT_REPOSITORY ${CINT_SRC}
    GIT_TAG ${CINT_VER}
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    PREFIX ${PROJECT_BINARY_DIR}/deps
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=RELEASE
        -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/deps
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DWITH_RANGE_COULOMB:STRING=1
        -DWITH_FORTRAN:STRING=0
        -DWITH_CINT2_INTERFACE:STRING=1
        -DWITH_F12:STRING=1
        -DWITH_4C1E:STRING=1
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/deps/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
