cmake_minimum_required(VERSION 3.13)
project(GRu_Mo)

set(CMAKE_CXX_STANDARD 17)
#TODO: replace uncomment version with version that will apply optimisations
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -l sqlite3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -l sqlite3")
add_executable(GRu_Mo main.cpp)