// a0 = 14.095
// a1 = -7.6842
// a2 = 2.8190
// a3 = -1.5369
// a4 = 0.82913
// a5 = -0.59109

#include <iostream>
#include "chebyshev_polynomials.hpp"
#include "../common/consts.hpp"

long double find_pol_in_x(const std::vector<long long> &T, long double x) {
    long double ans = 0;
    long double x_i = 1;
    for (long long coef: T) {
        ans += x_i * coef;
        x_i *= x;
    }
    return ans;
}

std::vector<long double> solve(unsigned N) {
    std::vector<long double> ans(N + 1); // k \in [0, N]
    adaai::ChebyshevPolynomials pols(N + 1);
    for (unsigned k = 0; k <= N; k++) {
        for (unsigned i = 1; i <= N + 1; i++) {
            long double x_i = cos((((double) (i * 2 - 1)) / ((double) (2 * N + 2)) * 3.14159265358979323846));
            ans[k] = std::exp(acos(x_i)) * find_pol_in_x(pols.get_polynomial(k), x_i);
            ans[k] = (ans[k] * 2) / (N + 1);
        }
    }
    return ans;
}

int main() {
    std::vector<long double> a = solve(20);
    std::cout.precision(17);
    for(long double a_k: a){
        std::cout << a_k << ",\n";
    }
    return 0;
}