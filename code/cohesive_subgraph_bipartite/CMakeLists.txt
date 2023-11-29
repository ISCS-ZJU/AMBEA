cmake_minimum_required(VERSION 3.10)
project(mbbp)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")

#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-march=native")
#find_package( Boost COMPONENTS filesystem system REQUIRED )

#include_directories(
#        ${BOOST_INCLUDE_DIRS}
#)
#add_compile_options(-march=native)

set(CMAKE_CXX_STANDARD 14)





add_executable(mbbp main.cpp bgraph.cpp Utility.cpp Utility.h bgraph.h basicEnu.cpp basicEnu.h denseMBB.cpp denseMBB.h LCUtility.cpp LCUtility.h basicEnuSet.cpp basicEnuSet.h UniCore.cpp UniCore.h BiCore.cpp BiCore.h mBBE.cpp mBBE.h ooMBE.cpp ooMBE.h)
