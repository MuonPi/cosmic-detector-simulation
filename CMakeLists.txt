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
#ADD_DEFINITIONS(-ffast-math)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(SOURCE_FILES
    "main.cpp"
    "algebra_utils.cpp"
    "geometry_types.cpp"
    "histogram.cpp"
)

set(HEADER_FILES
    "histogram.h"
    "sampled_distribution.h"
    "algebra_types.h"
    "algebra_utils.h"
    "geometry_types.h"
    "detector_setup.h"
)

# tell cmake to build our executable
ADD_EXECUTABLE(detsim 
    ${SOURCE_FILES}
    ${HEADER_FILES}
)

set(ALL_FILES
    "${SOURCE_FILES}"
    "${HEADER_FILES}"
)

TARGET_LINK_LIBRARIES(detsim
 pthread
)

add_custom_target(clangformat COMMAND clang-format-11 -style=WebKit -i ${ALL_FILES})
