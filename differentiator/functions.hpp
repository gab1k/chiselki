#ifndef F1_HPP_
#define F1_HPP_

#include "enums.hpp"
#include "../common/consts.hpp"

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

struct Func3 {
    template<typename T>
    T operator()(T x, T y) const { // Func3(x, y) = (8x^2 - 5xy - 3)/(x^2 + y^4 + 7)

        T xy5 = x;
        x *= y;
        xy5 *= 5;

        T x2_8 = x;
        x2_8 *= x;
        x2_8 *= 8;

        x2_8 -= xy5;
        x2_8 -= 3; // числитель

        T x2 = x * x;

        T y4 = y * y * y * y;

        x2 += y4;
        x2 += 7; // знаминатель

        x2_8 /= x2;

        return x2_8;
    }

    static double get_diff(WhichD T, double x, double y) {
        double res;
        if (T == WhichD::x) {
            res = (5 * pow(x, 2) * y + 2 * x * (8 * pow(y, 4) + 59) - 5 * y * (pow(y, 4) + 7)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 2);
        } else if (T == WhichD::y) {
            res = (4 * pow(y, 3) * (-8 * pow(x, 2) + 5 * x * y + 3) - 5 * x * (pow(x, 2) + pow(y, 4) + 7)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 2);
        } else if (T == WhichD::xx) {
            res = (-10 * pow(x, 3) * y - 6 * pow(x, 2) * (8 * pow(y, 4) + 59) + 30 * x * (pow(y, 4) + 7) * y +
                   16 * pow(y, 8) + 230 * pow(y, 4) + 826) / pow(pow(x, 2) + pow(y, 4) + 7, 3);
        } else if (T == WhichD::yy) {
            res = -(4 * pow(y, 2) * (24 * pow(x, 4) - 25 * pow(x, 3) * y + pow(x, 2) * (159 - 40 * pow(y, 4)) +
                                     5 * x * y * (3 * pow(y, 4) - 35) + 15 * pow(y, 4) - 63)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 3);
        } else {
            res = (5 * pow(x, 4) + 64 * pow(x, 3) * pow(y, 3) - 60 * pow(x, 2) * pow(y, 4) -
                   16 * x * pow(y, 3) * (4 * pow(y, 4) + 31) + 5 * (3 * pow(y, 8) + 14 * pow(y, 4) - 49)) /
                  pow(pow(x, 2) + pow(y, 4) + 7, 3);
        }
        return res;
    }
};

#endif // F1_HPP_
