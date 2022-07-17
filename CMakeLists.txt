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

set(SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/algebra_utils.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/geometry_types.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/histogram.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/detector_setup.cpp"
)

set(HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/histogram.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/sampled_distribution.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/algebra_types.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/algebra_utils.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/geometry_types.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/detector_setup.h"
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
