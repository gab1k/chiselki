#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>

#include "aad22.hpp"
#include "consts.hpp"
#include "differentiator.hpp"
#include "random_num.hpp"


std::unordered_map<WhichD, std::string> namesT{
        {WhichD::x,  "X"},
        {WhichD::y,  "Y"},
        {WhichD::xx, "XX"},
        {WhichD::yy, "YY"},
        {WhichD::xy, "XY"}
};
std::unordered_map<DiffMethod, std::string> namesM{
        {DiffMethod::stencil3,      "stencil3"},
        {DiffMethod::stencil3Extra, "stencil3Extra"},
        {DiffMethod::stencil5,      "stencil5"},
        {DiffMethod::stencil5Extra, "stencil5Extra"},
        {DiffMethod::FwdAAD,        "FwdAAD"}
};

struct Func1 {
    template<typename T>
    T operator()(T x, T y) const {
        return sin(x) * (exp(x * y) + cos(3 * x)) / (y * y + 1) - (adaai::Pi<double> * sin(y)) / (2 * exp(x));
    }

    static double get_diff(WhichD T, double x, double y) {
        double res, pi = adaai::Pi<double>;
        if (T == WhichD::x) {
            res = pi / (2 * exp(x)) * sin(y) + (y * exp(x * y) - 3 * sin(3 * x)) * sin(x) / (y * y + 1) +
                  (exp(x * y) + cos(3 * x)) * cos(x) / (y * y + 1);
        } else if (T == WhichD::y) {
            res = x * exp(x * y) * sin(x) / (y * y + 1) -
                  2 * y * (exp(x * y) + cos(3 * x)) * sin(x) / pow((y * y + 1), 2) -
                  pi * cos(y) / (2 * exp(x));
        } else if (T == WhichD::xx) {
            res = -pi * sin(y) / (2 * exp(x)) + 2 * (y * exp(x * y) - 3 * sin(3 * x)) * cos(x) / (y * y + 1) +
                  (y * y * exp(x * y) - 9 * cos(3 * x)) * sin(x) / (y * y + 1) -
                  (exp(x * y) + cos(3 * x)) * sin(x) / (y * y + 1);
        } else if (T == WhichD::yy) {
            res = x * x * exp(x * y) * sin(x) / (y * y + 1) - 4 * x * y * exp(x * y) * sin(x) / pow(y * y + 1, 2) +
                  8 * y * y * (exp(x * y) + cos(3 * x)) * sin(x) / pow(y * y + 1, 3) + pi * sin(y) / (2 * exp(x)) -
                  2 * (exp(x * y) + cos(3 * x)) * sin(x) / pow(y * y + 1, 2);
        } else {
            res = x * exp(x * y) * cos(x) / (y * y + 1) -
                  2 * y * (y * exp(x * y) - 3 * sin(3 * x)) * sin(x) / pow(y * y + 1, 2) -
                  2 * y * (exp(x * y) + cos(3 * x)) * cos(x) / pow(y * y + 1, 2) + pi * cos(y) / (2 * exp(x)) +
                  (x * y * exp(x * y) + exp(x * y)) * sin(x) / (y * y + 1);
        }
        return res;
    }
};

struct Func2 {
    template<typename T>
    T operator()(T x, T y) const {

        T x3 = x;
        x3 *= 3;

        T yy = y;
        yy *= y;
        yy += 1;

        T xy = x;
        xy *= y;

        return sin(x) * (exp(xy) + cos(x3)) / (yy) - (adaai::Pi<double> * sin(y)) / (2 * exp(x));
    }

    static double get_diff(WhichD T, double x, double y) {
        double res, pi = adaai::Pi<double>;
        if (T == WhichD::x) {
            res = pi / (2 * exp(x)) * sin(y) +
                  (y * exp(x * y) - 3 * sin(3 * x)) * sin(x) / (y * y + 1) +
                  (exp(x * y) + cos(3 * x)) * cos(x) / (y * y + 1);
        } else if (T == WhichD::y) {
            res = x * exp(x * y) * sin(x) / (y * y + 1) -
                  2 * y * (exp(x * y) + cos(3 * x)) * sin(x) / pow((y * y + 1), 2) -
                  pi * cos(y) / (2 * exp(x));
        } else if (T == WhichD::xx) {
            res = -pi * sin(y) / (2 * exp(x)) +
                  2 * (y * exp(x * y) - 3 * sin(3 * x)) * cos(x) / (y * y + 1) +
                  (y * y * exp(x * y) - 9 * cos(3 * x)) * sin(x) / (y * y + 1) -
                  (exp(x * y) + cos(3 * x)) * sin(x) / (y * y + 1);
        } else if (T == WhichD::yy) {
            res = x * x * exp(x * y) * sin(x) / (y * y + 1) -
                  4 * x * y * exp(x * y) * sin(x) / pow(y * y + 1, 2) +
                  8 * y * y * (exp(x * y) + cos(3 * x)) * sin(x) / pow(y * y + 1, 3) +
                  pi * sin(y) / (2 * exp(x)) -
                  2 * (exp(x * y) + cos(3 * x)) * sin(x) / pow(y * y + 1, 2);
        } else {
            res = x * exp(x * y) * cos(x) / (y * y + 1) -
                  2 * y * (y * exp(x * y) - 3 * sin(3 * x)) * sin(x) / pow(y * y + 1, 2) -
                  2 * y * (exp(x * y) + cos(3 * x)) * cos(x) / pow(y * y + 1, 2) +
                  pi * cos(y) / (2 * exp(x)) +
                  (x * y * exp(x * y) + exp(x * y)) * sin(x) / (y * y + 1);
        }
        return res;
    }
};

struct Func3 {
    template<typename T>
    T operator()(T x, T y) const { // Func3(x, y) = (8x^2 - 5xy - 3)/(x^2 + y^4 + 7)
        T x2_8 = x;
        x2_8 *= x;
        x2_8 *= 8;

        T xy5 = x;
        xy5 *= y;
        xy5 *= 5;


        x2_8 -= xy5;
        x2_8 -= 3; // числитель

        T x2 = x * x;

        T y4 = y * y * y * y;

        x2 += y4;
        x2 += 7; // знаменатель

        x2_8 /= x2;

        return x2_8;
//        return (8 * x * x - 5 * x * y - 3) / (x * x + y * y * y * y + 7);
    }

    static double get_diff(WhichD T, double x, double y) {
        double res;
        if (T == WhichD::x) {
            res = (5 * pow(x, 2) * y + 2 * x * (8 * pow(y, 4) + 59) - 5 * y * (pow(y, 4) + 7)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 2);
        } else if (T == WhichD::y) {
            res = (4 * pow(y, 3) * (-8 * pow(x, 2) + 5 * x * y + 3) -
                   5 * x * (pow(x, 2) + pow(y, 4) + 7)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 2);
        } else if (T == WhichD::xx) {
            res = (-10 * pow(x, 3) * y - 6 * pow(x, 2) * (8 * pow(y, 4) + 59) +
                   30 * x * (pow(y, 4) + 7) * y +
                   16 * pow(y, 8) + 230 * pow(y, 4) + 826) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 3);
        } else if (T == WhichD::yy) {
            res = -(4 * pow(y, 2) *
                    (24 * pow(x, 4) - 25 * pow(x, 3) * y + pow(x, 2) * (159 - 40 * pow(y, 4)) +
                     5 * x * y * (3 * pow(y, 4) - 35) + 15 * pow(y, 4) - 63)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 3);
        } else {
            res = (5 * pow(x, 4) + 64 * pow(x, 3) * pow(y, 3) -
                   60 * pow(x, 2) * pow(y, 4) -
                   16 * x * pow(y, 3) * (4 * pow(y, 4) + 31) +
                   5 * (3 * pow(y, 8) + 14 * pow(y, 4) - 49)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 3);
        }
        return res;
    }
};


template<typename Callable>
double get_res(Callable F, WhichD T, DiffMethod M, double x, double y) {
    if (T == WhichD::x) {
        if (M == DiffMethod::stencil3) {
            return Differentiator<WhichD::x, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return Differentiator<WhichD::x, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return Differentiator<WhichD::x, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return Differentiator<WhichD::x, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return Differentiator<WhichD::x, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::y) {
        if (M == DiffMethod::stencil3) {
            return Differentiator<WhichD::y, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return Differentiator<WhichD::y, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return Differentiator<WhichD::y, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return Differentiator<WhichD::y, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return Differentiator<WhichD::y, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::xx) {
        if (M == DiffMethod::stencil3) {
            return Differentiator<WhichD::xx, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return Differentiator<WhichD::xx, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return Differentiator<WhichD::xx, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return Differentiator<WhichD::xx, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return Differentiator<WhichD::xx, DiffMethod::FwdAAD>(F, x, y);
        }
    } else if (T == WhichD::yy) {
        if (M == DiffMethod::stencil3) {
            return Differentiator<WhichD::yy, DiffMethod::stencil3>(F, x, y);
        } else if (M == DiffMethod::stencil5) {
            return Differentiator<WhichD::yy, DiffMethod::stencil5>(F, x, y);
        } else if (M == DiffMethod::stencil3Extra) {
            return Differentiator<WhichD::yy, DiffMethod::stencil3Extra>(F, x, y);
        } else if (M == DiffMethod::stencil5Extra) {
            return Differentiator<WhichD::yy, DiffMethod::stencil5Extra>(F, x, y);
        } else if (M == DiffMethod::FwdAAD) {
            return Differentiator<WhichD::yy, DiffMethod::FwdAAD>(F, x, y);
        }
    }
    if (M == DiffMethod::stencil3) {
        return Differentiator<WhichD::xy, DiffMethod::stencil3>(F, x, y);
    } else if (M == DiffMethod::stencil5) {
        return Differentiator<WhichD::xy, DiffMethod::stencil5>(F, x, y);
    } else if (M == DiffMethod::stencil3Extra) {
        return Differentiator<WhichD::xy, DiffMethod::stencil3Extra>(F, x, y);
    } else if (M == DiffMethod::stencil5Extra) {
        return Differentiator<WhichD::xy, DiffMethod::stencil5Extra>(F, x, y);
    }
    return Differentiator<WhichD::xy, DiffMethod::FwdAAD>(F, x, y);
}

template<typename Callable>
std::map<std::pair<DiffMethod, WhichD>, double>
stress_test(Callable F, double from_x, double to_x, double from_y, double to_y, int n = 10000) {
    std::map<std::pair<DiffMethod, WhichD>, double> res;
    for (const auto &TT: namesT) {
        for (const auto &MM: namesM) {
            res[{MM.first, TT.first}] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        double x = get_rand_val(from_x, to_x);
        double y = get_rand_val(from_y, to_y);
        for (const auto &TT: namesT) {
            for (const auto &MM: namesM) {
                double expected = F.get_diff(TT.first, x, y);
                double current = get_res(F, TT.first, MM.first, x, y);
                res[{MM.first, TT.first}] = std::max(res[{MM.first, TT.first}], std::abs(expected - current));
            }
        }
    }
    return res;
}


int main() {
    std::string filePath = "../differentiator/logs/logs.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out);

    Func1 f1;
    ofs << "Testing Func1:\n";
    auto res = stress_test(f1, -3, 3, -3, 3);
    for (const auto &MM: namesM) {
        ofs << "Diff method " << MM.second << ":\n";
        for (const auto &TT: namesT) {
            ofs << "Diff type: " << TT.second << " max error is: " << res[{MM.first, TT.first}] << "\n";
        }
        ofs << "\n";
    }

    ofs << "\n\n";

    Func2 f2; // Func2 == Func1, but with +=, *= etc
    ofs << "Testing Func2:\n";
    res = stress_test(f2, -3, 3, -3, 3);
    for (const auto &MM: namesM) {
        ofs << "Diff method " << MM.second << ":\n";
        for (const auto &TT: namesT) {
            ofs << "Diff type: " << TT.second << " max error is: " << res[{MM.first, TT.first}] << "\n";
        }
        ofs << "\n";
    }
    ofs << "\n\n";

    Func3 f3;
    ofs << "Testing Func3:\n";
    res = stress_test(f3, -1000, 1000, -1000, 1000, 10000);
    for (const auto &MM: namesM) {
        ofs << "Diff method " << MM.second << ":\n";
        for (const auto &TT: namesT) {
            ofs << "Diff type: " << TT.second << " max error is: " << res[{MM.first, TT.first}] << "\n";
        }
        ofs << "\n";
    }
    ofs << "\n\n";
    ofs.close();
    return 0;
}
