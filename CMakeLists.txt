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
if(NOT memalloc_FOUND)
    ExternalProject_Add(
        Memalloc
        SOURCE_DIR /sdcard/Jay/Projects/xab_laud_ubject/memalloc/
        INSTALL_DIR ${CMAKE_BINARY_DIR}/memalloc
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        BUILD_ALWAYS OFF
        STAMP_DIR ${CMAKE_BINARY_DIR}/stamps
    )
endif()

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/zot/lib/cmake/zot")
find_package(zot QUIET)
if(NOT zot_FOUND)
    ExternalProject_Add(
        Zot
        SOURCE_DIR /sdcard/Jay/Projects/xab_laud_ubject/zot/
        INSTALL_DIR ${CMAKE_BINARY_DIR}/zot
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
        BUILD_ALWAYS OFF
        STAMP_DIR ${CMAKE_BINARY_DIR}/stamps
    )
    if(NOT memalloc_FOUND)
        add_dependencies(Zot Memalloc)
    endif()
endif()

add_library(lray STATIC)
include(sources.cmake)
include(headers.cmake)
target_sources(lray PRIVATE ${SOURCE_FILES})
if(NOT zot_FOUND)
    add_dependencies(lray Zot)
endif()
target_include_directories(lray PUBLIC src src/renderer ${CMAKE_BINARY_DIR}/zot/include)

add_executable(lray_exe main.c)
target_link_libraries(lray_exe PUBLIC memalloc zot lray)
if(NOT memalloc_FOUND)
    target_link_directories(lray_exe PRIVATE ${CMAKE_BINARY_DIR}/memalloc/lib)
    add_dependencies(lray_exe Memalloc)
endif()
if(NOT zot_FOUND)
    target_link_directories(lray_exe PRIVATE ${CMAKE_BINARY_DIR}/zot/lib)
    add_dependencies(lray_exe Zot)
endif()
add_dependencies(lray_exe lray)


install(TARGETS lray DESTINATION lib)
install(TARGETS lray_exe DESTINATION bin)


