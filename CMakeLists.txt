cmake_minimum_required(VERSION 3.27)
project(chiselki)

set(CMAKE_CXX_STANDARD 20)

include_directories(common)

include_directories(exp)
add_executable(test-exp
        common/consts.hpp
        common/doctest.h
        exp/exp.hpp
        exp/test_exp.cpp)

include_directories(useful)
add_executable(pre-calculation-chebyshev-coefficients
            useful/chebyshev_coefficients.cpp)

find_package(GSL REQUIRED)
target_link_libraries(pre-calculation-chebyshev-coefficients GSL::gsl)
