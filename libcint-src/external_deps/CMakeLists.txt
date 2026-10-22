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

# CINT_SRC is cloned with git, unless it is a local path that is not a git
# repository (a source archive or an unpacked source directory), or a URL of a
# source archive. CINT_VER applies to git sources only, CINT_SRC_SHA256 to
# archives only.
if((EXISTS "${CINT_SRC}" AND NOT EXISTS "${CINT_SRC}/.git")
    OR CINT_SRC MATCHES "\\.(tar\\.gz|tgz|tar\\.bz2|tar\\.xz|zip)$")
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
    endif()
    set(CINT_DOWNLOAD URL ${CINT_SRC})
    if(CINT_SRC_SHA256)
        list(APPEND CINT_DOWNLOAD URL_HASH SHA256=${CINT_SRC_SHA256})
    endif()
else()
    if(CINT_SRC_SHA256)
        message(FATAL_ERROR "CINT_SRC_SHA256 applies to archive sources only, not git: ${CINT_SRC}")
    endif()
    set(CINT_DOWNLOAD GIT_REPOSITORY ${CINT_SRC} GIT_TAG ${CINT_VER} GIT_SHALLOW TRUE GIT_PROGRESS TRUE)
endif()

ExternalProject_Add(cint
    ${CINT_DOWNLOAD}
    PREFIX ${PROJECT_BINARY_DIR}/deps
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=RELEASE
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/deps
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DWITH_RANGE_COULOMB:STRING=1
        -DWITH_FORTRAN:STRING=0
        -DWITH_CINT2_INTERFACE:STRING=1
        -DWITH_F12:STRING=${WITH_F12}
        -DWITH_4C1E:STRING=${WITH_4C1E}
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/deps/${CMAKE_INSTALL_LIBDIR} DESTINATION ${CMAKE_INSTALL_PREFIX})
