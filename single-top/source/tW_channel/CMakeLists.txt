cmake_minimum_required(VERSION 3.10)
project(withtW_analysis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

# ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# Executável
add_executable(run_withtW
    top_inclusive_withBH_withtW.cpp)

target_link_libraries(run_withtW ${ROOT_LIBRARIES})

