cmake_minimum_required(VERSION 3.10.1)
project(MBEA_proj_interface)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-L/usr/lib -Wl,--no-as-needed,-lprofiler,--as-needed")

set(SOURCE_FILES
        BicliqueFinder.cc
        BicliqueFinder.h
        BiGraph.cc
        BiGraph.h 
        utility.cc
        utility.h main.cc)

add_executable(MBEA_proj_interface ${SOURCE_FILES})
