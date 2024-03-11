#include "differentiator.hpp"
#include "aad22.cpp"
#include <iostream>
#include <vector>
#include <unordered_map>
#include "F1.hpp"
#include "../common/random_num.hpp"

template<typename Callable>
double get_res(Callable F, WhichD T, DiffMethod M, double x, double y) {
    if (T == WhichD::x) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return adaai::Differentiator<WhichD::x, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::y) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return adaai::Differentiator<WhichD::y, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::xx) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::yy) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::FwdAAD>(F, x, y);
        }
    }
    if (M == DiffMethod::stencil3) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil3>(F, x, y);
    } else if (M == DiffMethod::stencil5) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil5>(F, x, y);
    } else if (M == DiffMethod::stencil3Extra) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil3Extra>(F, x, y);
    } else if (M == DiffMethod::stencil5Extra) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil5Extra>(F, x, y);
    }
    return adaai::Differentiator<WhichD::xy, DiffMethod::FwdAAD>(F, x, y);
}

template<typename Callable>
std::unordered_map<DiffMethod, double>
stress_test(Callable F, double from_x, double to_x, double from_y, double to_y, int n = 10000) {
    std::unordered_map<DiffMethod, double> res{{DiffMethod::stencil3,      0},
                                               {DiffMethod::stencil5,      0},
                                               {DiffMethod::stencil3Extra, 0},
                                               {DiffMethod::stencil5Extra, 0},
                                               {DiffMethod::FwdAAD,        0}};
    for (int i = 0; i < n; i++) {
        double x = get_rand_val(from_x, to_x);
        double y = get_rand_val(from_y, to_y);
        for (const auto &TT: namesT) {
            for (const auto &MM: namesM) {
                double expected = F.get_diff(TT.first, x, y);
                double current = get_res(F, TT.first, MM.first, x, y);
                res[MM.first] = std::max(res[MM.first], std::abs(expected - current));
            }
        }
    }
    return res;
}


int main() {
    Func1 f;
    auto res = stress_test(f, -3, 3, -3, 3);
    for(auto P: res){
        std::cout << "Method " << namesM[P.first] << " max error is: " << P.second << "\n";
    }
    return 0;
}