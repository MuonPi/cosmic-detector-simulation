cmake_minimum_required(VERSION 2.6)

PROJECT(DetectorSimulator)

include(GNUInstallDirs)

#SET(CMAKE_VERBOSE_MAKEFILE ON)

#The following lines are required to use Dart
#ENABLE_TESTING()
#INCLUDE(Dart)

SET(CMAKE_BUILD_TYPE Debug)

ADD_DEFINITIONS(-std=gnu++17)
ADD_DEFINITIONS(-O3)
ADD_DEFINITIONS(-ffast-math)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# tell cmake to build our executable
ADD_EXECUTABLE(detsim 
    algebra_utils.cpp
    geometry_types.cpp
    histogram.cpp
    main.cpp
)

TARGET_LINK_LIBRARIES(detsim
 pthread
)
