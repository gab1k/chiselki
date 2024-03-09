#include "differentiator.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>

double F1(double x, double y) {
    return exp(3 * x) / 13 + cos(x * y) * sin(x);
}

// dF1/dx = -y * sin(x) * sin(x*y) + 3/13 * e^{3x} + cos(x) * cos(xy)
// dF1/dy = -x * sin(x) * sin(xy)
// d^2 F1/dxdx = -y^2 * sin(x) * cos(xy) - 2y * sin(xy) * cos(x) + 9/13 * e^{3x} - sin(x) * cos(xy)
// d^2 F1/dydy = -x^2 * sin(x) * cos(xy)
// d^2 F1/dxdy = xy * sin(x) * cos(xy) - x * sin(xy) * cos(x) - sin(x) * sin(xy)

AAD22 FADD22(AAD22 x, AAD22 y){
    return 7 * x * y * x;
}


double get_resF1(WhichD T, DiffMethod M, double x, double y) {
    if (T == WhichD::x) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil3>(F1, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil5>(F1, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil3Extra>(F1, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::x, DiffMethod::stencil5Extra>(F1, x, y);
        }
    } else if (T == WhichD::y) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil3>(F1, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil5>(F1, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil3Extra>(F1, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::y, DiffMethod::stencil5Extra>(F1, x, y);
        }
    } else if (T == WhichD::xx) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil3>(F1, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil5>(F1, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil3Extra>(F1, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::xx, DiffMethod::stencil5Extra>(F1, x, y);
        }
    } else if (T == WhichD::yy) {
        if (M == DiffMethod::stencil3) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil3>(F1, x, y);
        } else if (M == DiffMethod::stencil5) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil5>(F1, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil3Extra>(F1, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return adaai::Differentiator<WhichD::yy, DiffMethod::stencil5Extra>(F1, x, y);
        }
    }
    if (M == DiffMethod::stencil3) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil3>(F1, x, y);
    } else if (M == DiffMethod::stencil5) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil5>(F1, x, y);
    } else if (M == DiffMethod::stencil3Extra) {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil3Extra>(F1, x, y);
    } else {
        return adaai::Differentiator<WhichD::xy, DiffMethod::stencil5Extra>(F1, x, y);
    }
}

double F2(double x, double y) {
    return (y + 2) * (sin(x) - cos(y)) - (x * 7) / exp(y);
}

double F_easy(double x, double y) {
    return 7 * x * y;
}

int main() {

    std::vector<std::pair<double, double>> points{{0, 0},
                                                  {1, 3},
    };
    std::unordered_map<WhichD, std::vector<double>> expect_f1 = {
            {WhichD::x,  {16.0 / 13.0,
                                 cos(1) * cos(3) - sin(1) * sin(3) + (3 * pow(exp(1), 3)) / 13}},
            {WhichD::y,  {0,
                                 -sin(1) * sin(3)}},
            {WhichD::xx, {9.0 / 13.0,
                                 -6 * sin(3) * cos(1) - 10 * sin(1) * cos(3) + (9 * pow(exp(1), 3)) / 13}},
            {WhichD::yy, {0,
                                 -sin(1) * cos(3)}},
            {WhichD::xy, {0,
                                 -sin(1) * sin(3) - sin(3) * cos(1) - 3 * sin(1) * cos(3)}},
    };
    std::unordered_map<WhichD, std::string> namesT{
            {WhichD::x, "X"}, {WhichD::y, "Y"}, {WhichD::xx, "XX"}, {WhichD::yy, "YY"}, {WhichD::xy, "XY"}
    };
    std::unordered_map<DiffMethod, std::string> namesM{
            {DiffMethod::stencil3, "stencil3"}, {DiffMethod::stencil3Extra, "stencil3Extra"},
            {DiffMethod::stencil5, "stencil5"}, {DiffMethod::stencil5Extra, "stencil5Extra"},
            {DiffMethod::FwdAAD, "FwdAAD"}
    };

    std::vector<WhichD> diff_type{WhichD::x, WhichD::y, WhichD::xx, WhichD::yy, WhichD::xy};
    std::vector<DiffMethod> diff_method{DiffMethod::stencil3, DiffMethod::stencil3Extra, DiffMethod::stencil5,
                                        DiffMethod::stencil5Extra};

    for (WhichD T: diff_type) {
        for (DiffMethod M: diff_method) {
            for (int i = 0; i < points.size(); i++) {
                std::cout << "Diff type - " << namesT[T] << "\n" << "Diff Method - " << namesM[M] << "\n";
                std::cout << "Expected: " << expect_f1[T][i] << "\n";
                std::cout << "Result:   " << get_resF1(T, M, points[i].first, points[i].second) << "\n"
                          << "\n\n";

            }
        }
    }
    return 0;
}