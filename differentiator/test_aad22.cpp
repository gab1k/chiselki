#include "../common/random_num.hpp"
#include "differentiator.hpp"
#include "aad22.cpp"
#include "F1.hpp"

#include <iostream>

template<typename Callable>
double get_res_fwd(Callable F, WhichD T, double x, double y) {
    if (T == WhichD::x) {
        return adaai::Differentiator<WhichD::x, DiffMethod::FwdAAD>(F, x, y);
    } else if (T == WhichD::y) {
        return adaai::Differentiator<WhichD::y, DiffMethod::FwdAAD>(F, x, y);
    } else if (T == WhichD::xx) {
        return adaai::Differentiator<WhichD::xx, DiffMethod::FwdAAD>(F, x, y);
    } else if (T == WhichD::yy) {
        return adaai::Differentiator<WhichD::yy, DiffMethod::FwdAAD>(F, x, y);
    }
    return adaai::Differentiator<WhichD::xy, DiffMethod::FwdAAD>(F, x, y);
}

template<typename Callable>
double stress_test(Callable F, double from_x, double to_x, double from_y, double to_y, int n = 10000) {
    double max_err = 0;
    for (int i = 0; i < n; i++) {
        double x = get_rand_val(from_x, to_x);
        double y = get_rand_val(from_y, to_y);
//        std::cout << "Point is: {" << x << ", " << y << "}\n";
        for (const auto &TT: namesT) {
            double expected = F.get_diff(TT.first, x, y);
            double current = get_res_fwd(F, TT.first, x, y);
//            std::cout << "expected: " << expected << "\n";
//            std::cout << "current:  " << current << "\n\n";
            max_err = std::max(max_err, std::abs(expected - current));
        }
    }
    return max_err;
}


int main() {
    Func1 f;
    double max_err = stress_test(f, -5, 5, -5, 5);
    std::cout << "Max absolute error is: " << max_err << "\n";
    return 0;
}