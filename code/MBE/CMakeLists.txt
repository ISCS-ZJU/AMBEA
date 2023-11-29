cmake_minimum_required(VERSION 3.3.2)
project(MBE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-L/usr/lib -Wl,--no-as-needed,-lprofiler,--as-needed")

set(SOURCE_FILES
        BicliqueFinder.cpp
        BicliqueFinder.h
        BiGraph.cpp
        BiGraph.h 
        Utility.cpp
        Utility.h
        BaselineFinder.cpp
        BaselineFinder.h
        AggressiveFinder.cpp
        AggressiveFinder.h
        main.cpp)

set(TBB_ROOT "oneTBB")
include_directories(${TBB_ROOT}/include)
link_directories(${TBB_ROOT}/lib)

add_executable(MBE_ALL ${SOURCE_FILES})
target_link_libraries(MBE_ALL tbb)

# find_package(TBB REQUIRED)
# add_executable(MBE_ALL ${SOURCE_FILES})
# target_link_libraries(MBE_ALL TBB::tbb)
