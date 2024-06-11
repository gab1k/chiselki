cmake_minimum_required(VERSION 3.27)
project(chiselki)

set(CMAKE_CXX_STANDARD 20)

find_package(GSL REQUIRED)


include_directories(useful)


include_directories(exp)
add_executable(test-exp
        useful/consts.hpp
        useful/doctest.h
        exp/exp.hpp
        exp/exp_taylor.hpp
        exp/exp_pade.hpp
        exp/exp_chebyshev.hpp
        exp/test_exp.cpp)
target_link_libraries(test-exp GSL::gsl)

add_executable(get-chebyshev-coefficients
        exp/chebyshev_coefficients.cpp)
target_link_libraries(get-chebyshev-coefficients GSL::gsl)


include_directories(exp/fft)
add_executable(test-exp-fft
        useful/consts.hpp
        useful/doctest.h
        exp/fft/exp_fft.hpp
        exp/fft/test_exp_fft.cpp)
target_link_libraries(test-exp-fft GSL::gsl)

add_executable(get-fft-coefficients
        exp/fft/chebyshev_polynomials.hpp
        exp/fft/fft_coefficients.cpp)
target_link_libraries(get-fft-coefficients GSL::gsl)


include_directories(differentiator)
add_executable(test-differentiator
        useful/consts.hpp
        useful/random_num.hpp
        differentiator/differentiator.hpp
        differentiator/aad22.hpp
        differentiator/aad22.cpp
        differentiator/test_differentiator.cpp)


include_directories(atmosphere_model)
add_executable(ode
        atmosphere_model/rho.hpp
        atmosphere_model/rho.cpp
        atmosphere_model/ode_solver.cpp)

add_executable(get-rho-coefficients
        atmosphere_model/rho.hpp
        atmosphere_model/rho.cpp
        atmosphere_model/rho_coefficients.cpp)


#include_directories(everhart_method)
#add_executable(evethart
#        everhart_method/stepper.hpp
#        everhart_method/determine_satellite_position.cpp)
