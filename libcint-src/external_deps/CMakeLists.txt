# Please note that options (`CINT_SRC`, `CINT_VER`, etc.) will be always be overrided
# in build.rs building script. However, this CMakeLists.txt is also valid for standalone
# building.

cmake_minimum_required(VERSION 3.15)
project(dep-libcint)

include(ExternalProject)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

if(CINT_SRC)
    message(STATUS "User specified directory of source code (libcint): ${CINT_SRC}")
else()
    set(CINT_SRC https://github.com/sunqm/libcint.git)
    message(STATUS "Download source code (libcint): ${CINT_SRC}")
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
        -DWITH_F12:STRING=${WITH_F12}
        -DWITH_4C1E:STRING=${WITH_4C1E}
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/deps/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
