cmake_minimum_required(VERSION 2.8)
project(geohash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

enable_testing()

add_executable(test_geohash geohash.cpp test_geohash.cpp)

add_test(geohash test_geohash)
