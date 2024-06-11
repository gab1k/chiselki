#include <cmath>
#include <iostream>
#include <vector>

#include "rho.hpp"

double get_Q(double y, double v) { // y = h
    double diam = 0.216;
    double S = (M_PI * diam * diam) / 4.0;
    Rho rho;
    CD cd;
    double M = v / (sqrt(rho.get_p(y) / rho(y)));
    return (cd(M) * rho(y) * v * v * S) / 2;
}

std::vector<double> get_diff(const std::vector<double> &u) {
    // 0: x, 1: v_x, 2: y, 3: v_y
    double v = sqrt(u[1] * u[1] + u[3] * u[3]);
    std::vector<double> res{u[1], -1, u[3], -1};
    const double g = 9.80655;
    double m = 106;
    double Q = get_Q(u[2], v);
    res[1] = -((Q * u[1]) / v) / m;
    res[3] = -((Q * u[3]) / v) / m - g;
    return res;
}

void find_a_b() {
    double p1 = 1.5, p2 = 2.0;
    // CD(M) = a * M^{-b}

    // CD(p1) = a * p1{-b}
    // CD(p2) = a * p2^{-b}
    // CD(p1) / CD(p2) = (p1/p2) ^ {-b} = (p2/p1)^b
    // log (CD(p1) / CD(p2)) = log ((p2/p1)^b) = b * log (p2/p1)
    // b = log (CD(p1) / CD(p2)) / log (p2/p1)
    CD cd;
    double cd1 = cd(p1);
    double cd2 = cd(p2);
    double b = log(cd1 / cd2) / log(p2 / p1);
    // a = CD(p1) / (p1^{-b})
    double a = cd1 / pow(p1, -b);

    std::cout.precision(17);
    std::cout << "a = " << a << "\nb = " << b << "\n";
}
