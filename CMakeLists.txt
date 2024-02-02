cmake_minimum_required(VERSION 3.27)
project(chiselki)

set(CMAKE_CXX_STANDARD 20)

include_directories(Exp)

add_executable(test-exp
        Exp/Consts.hpp
        Exp/Exp.hpp
        Exp/TestExp.cpp)
