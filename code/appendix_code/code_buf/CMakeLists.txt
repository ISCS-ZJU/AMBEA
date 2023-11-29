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
        PMBE.cpp
        PMBE.h
        BicliqueFinderFast.cpp
        BicliqueFinderFast.h
        BitMBE.cpp
        BitMBE.h
        ComMBE.cpp
        ComMBE.h
        QuasiMBE.cpp
        QuasiMBE.h
        BoundMBE.cpp
        BoundMBE.h
        main.cpp)

find_package(TBB REQUIRED)
add_executable(MBE_ALL ${SOURCE_FILES})
target_link_libraries(MBE_ALL TBB::tbb)
