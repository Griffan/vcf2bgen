cmake_minimum_required(VERSION 2.8.4)
project(libVcf)

add_subdirectory(libStatGenLite)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)

add_library(Vcf ${SOURCE_FILES})
#set_target_properties(Vcf PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/libVcf")
target_include_directories(Vcf PRIVATE libStatGenLite)
add_dependencies(Vcf libStatGenLite)
