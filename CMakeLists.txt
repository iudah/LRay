cmake_minimum_required(VERSION 3.15)

project(lray)

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -mfpu=neon")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -static-libsan -mfpu=neon")
set(C_WARNINGS "-Wall -Wextra -Wshadow -Wformat=2 -Wunused")

include(ExternalProject)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/memalloc/lib/cmake/memalloc")

find_package(memalloc QUIET)
add_library(memalloc_imported STATIC IMPORTED)
if(NOT TARGET memalloc::memalloc)
  message(STATUS "memalloc not found, adding as ExternalProject")
  ExternalProject_Add(
    MemAlloc
    GIT_REPOSITORY "https://github.com/iudah/memalloc.git"
    GIT_TAG origin/master
    GIT_SHALLOW TRUE
    INSTALL_DIR ${CMAKE_BINARY_DIR}/memalloc
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DBUILD_EXECUTABLE=OFF
    BUILD_ALWAYS OFF
    STAMP_DIR ${CMAKE_BINARY_DIR}/stamps
    UPDATE_DISCONNECTED TRUE
    BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/memalloc/lib/${CMAKE_STATIC_LIBRARY_PREFIX}memalloc${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )

  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/memalloc/include")

  add_dependencies(memalloc_imported MemAlloc)
  add_library(memalloc::memalloc ALIAS memalloc_imported)
endif()
set_target_properties(memalloc_imported PROPERTIES
  IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/memalloc/lib/${CMAKE_STATIC_LIBRARY_PREFIX}memalloc${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/memalloc/include"
)

find_package(zot QUIET)
add_library(zot_imported STATIC IMPORTED)
if(NOT TARGET zot::zot)
  message(STATUS "zot not found, adding as ExternalProject")
  ExternalProject_Add(
    zot
    GIT_REPOSITORY "https://github.com/iudah/zot.git"
    GIT_TAG origin/main
    GIT_SHALLOW TRUE
    INSTALL_DIR ${CMAKE_BINARY_DIR}/zot
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DBUILD_EXECUTABLE=OFF
    BUILD_ALWAYS OFF
    STAMP_DIR ${CMAKE_BINARY_DIR}/stamps
    UPDATE_DISCONNECTED TRUE
    BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/zot/lib/${CMAKE_STATIC_LIBRARY_PREFIX}zot${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )

  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/zot/include")

  add_dependencies(zot_imported zot MemAlloc)
  add_library(zot::zot ALIAS zot_imported)
endif()
set_target_properties(zot_imported PROPERTIES
  IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/zot/lib/${CMAKE_STATIC_LIBRARY_PREFIX}zot${CMAKE_STATIC_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/zot/include"
)

add_library(lray STATIC)
include(sources.cmake)
include(headers.cmake)
target_sources(lray PRIVATE ${SOURCE_FILES})

target_include_directories(lray PUBLIC src/math src/renderer src/stb_image ${CMAKE_BINARY_DIR}/zot/include)

add_executable(lray_exe main.c)
target_link_libraries(lray_exe PUBLIC memalloc::memalloc zot::zot lray)


add_dependencies(lray_exe lray)


install(TARGETS lray DESTINATION lib)
install(TARGETS lray_exe DESTINATION bin)


