cmake_minimum_required(VERSION 2.8.4)
project(libbgen)

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)

file(GLOB HEADER_FILES
        "genfile/*.hpp"
        "genfile/bgen/*.hpp"
        )

add_library(bgen ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(bgen PRIVATE ../zstd-1.3.2/lib/)
add_dependencies(bgen libzstd_static)