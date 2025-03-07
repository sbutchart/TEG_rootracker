cmake_minimum_required(VERSION 3.9)

project(TEG_ROOTRACKER)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
ROOT_GENERATE_DICTIONARY(TEG_ROOTRACKER.cpp LINKDEF LinkDef.h)

add_executable(TEG_ROOTRACKER TEG_ROOTRACKER.cpp)
target_link_libraries(TEG_ROOTRACKER ${ROOT_LIBRARIES})
