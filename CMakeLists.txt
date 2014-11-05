cmake_minimum_required(VERSION 2.8)
project(geohash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-deprecated-register")
    set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -stdlib=libc++")
endif(APPLE)

enable_testing()

add_executable(test_geohash geohash.cpp test_geohash.cpp)

add_test(geohash test_geohash)
