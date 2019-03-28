cmake_minimum_required(VERSION 3.10)
project(GRu_Mo)

set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -l sqlite3")
add_executable(GRu_Mo src/main.cpp)
