#include <iostream>

#include "rho.hpp"

int main() {
    Rho rho;
    rho.print_all_coeff(); // correct only for correct layer pref
    find_a_b();
    double h = 0;
    while (h < 47000) {
        h += 100;
        std::cout << "h = " << h << ", Rho = " << rho(h) << "\n";
    }
}
