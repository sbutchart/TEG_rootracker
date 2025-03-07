cmake_minimum_required(VERSION 3.9)

project(TEG_ROOTRACKER)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT)
#find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
#find_program(LSB_RELEASE_EXEC lsb_release)
#find_package(TBB REQUIRED)

include(${ROOT_USE_FILE})
#string(CONCAT CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -Wall -Wextra -Werror -Wno-error=deprecated-copy")
#grabs the distribution
#execute_process(COMMAND ${LSB_RELEASE_EXEC} -is
  #OUTPUT_VARIABLE LSB_RELEASE_ID_SHORT
  #OUTPUT_STRIP_TRAILING_WHITESPACE
#)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
ROOT_GENERATE_DICTIONARY(TEG_ROOTRACKER.cpp LINKDEF LinkDef.h)

#add_library(TEG_ROOTRACKER_LIB SHARED TEG_ROOTRACKER.cpp)
#target_link_libraries(TEG_ROOTRACKER_LIB ${ROOT_LIBRARIES})

add_executable(TEG_ROOTRACKER TEG_ROOTRACKER.cpp)
target_link_libraries(TEG_ROOTRACKER ${ROOT_LIBRARIES})
