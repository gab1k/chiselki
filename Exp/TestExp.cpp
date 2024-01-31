#include "Exp.hpp"

int main(){
    float num = -3.0;
    std::cout.precision(20);
    std::cout << ADAAI::Exp(num) << " " << std::exp(num);
}