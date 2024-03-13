#include <iostream>
#include <unordered_map>
#include <fstream>
#include <map>

#include "differentiator.hpp"
#include "functions.hpp"
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
    std::string filePath = "../differentiator/logs.txt";
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