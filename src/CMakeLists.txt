cmake_minimum_required(VERSION 3.5)

project(lray)

add_executable(lray main.cpp)

install(TARGETS lray DESTINATION bin)

include(sources.cmake)

include(headers.cmake)
