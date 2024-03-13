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

include_directories(differentiator)
add_executable(test-differentiator
        common/random_num.hpp
        common/consts.hpp
        differentiator/differentiator.hpp
        differentiator/functions.hpp
        differentiator/enums.hpp
        differentiator/aad22.hpp
        differentiator/aad22.cpp
        differentiator/test_differentiator.cpp)

include_directories(useful)
add_executable(pre-calculation-chebyshev-coefficients
        useful/chebyshev_coefficients.cpp)

find_package(GSL REQUIRED)
target_link_libraries(pre-calculation-chebyshev-coefficients GSL::gsl)
target_link_libraries(test-exp GSL::gsl)
