cmake_minimum_required (VERSION 3.5.1)
project (MetaBAT)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message("Installing ${CMAKE_BUILD_TYPE} MetaBAT into ${CMAKE_INSTALL_PREFIX}")

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(cmake/zlib.cmake)
include(cmake/htslib.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_definitions(-D_XOPEN_SOURCE=700)

set(PRE_CONFIGURE_FILE "metabat_version.h.in")
set(POST_CONFIGURE_FILE "metabat_version.h")
include(cmake/git-watcher.cmake)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using AppleClang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mtune=native")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
endif()

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src)

if (NOT NO_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

install(PROGRAMS runMetaBat.sh merge_depths.pl aggregateBinDepths.pl aggregateContigOverlapsByBin.pl
        DESTINATION bin/
       )
INSTALL(CODE "execute_process( \
    COMMAND ${CMAKE_COMMAND} -E create_symlink \
    metabat2 \
    ${CMAKE_INSTALL_PREFIX}/bin/metabat   \
    )"
)


