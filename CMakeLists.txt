cmake_minimum_required(VERSION 3.0.0)
project(brightandcontrast VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(OpenCV REQUIRED)

aux_source_directory(./src DIR_SRCS)
include_directories(./inc)
add_executable(brightandcontrast ${DIR_SRCS})

target_link_libraries(brightandcontrast ${OpenCV_LIBS})

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)