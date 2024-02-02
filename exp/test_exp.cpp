#include "exp.hpp"
#include <iostream>
#include <complex>

int main(){
    float num = -3.0;
    std::cout.precision(20);
    std::cout << adaai::exp<float>(num) << " " << std::exp(num);
}
