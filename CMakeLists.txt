cmake_minimum_required(VERSION 2.8.4)
project(vcf2bgen)

ADD_DEFINITIONS(-DHAVE_ZLIB=1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0")

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
add_subdirectory(libVcf)
add_subdirectory(zstd-1.3.2/build/cmake/)
add_subdirectory(src)
add_executable(vcf2bgen main.cpp)
target_link_libraries(vcf2bgen Vcf bgen libzstd_static libStatGenLite z)

add_dependencies(vcf2bgen bgen Vcf)

