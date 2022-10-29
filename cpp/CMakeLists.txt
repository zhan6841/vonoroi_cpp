# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.4.1)

project(Voronoi_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -g -W -O3")

add_executable(voronoi_test 
    main.cpp
    VonoroiAdapt.cpp VonoroiAdapt.hpp
    Vonoroi.cpp Vonoroi.hpp
)