#ifndef F1_HPP_
#define F1_HPP_

#include "enums.hpp"
#include "../common/consts.hpp"

struct Func1 {
    template<typename T>
    T operator()(T x, T y) const{
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

#endif // F1_HPP_
