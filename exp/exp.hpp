#ifndef EXP_HPP_
#define EXP_HPP_

#include "consts.hpp"
#include <cmath>
#include <cstdint>

enum class MethodE : int {
    Taylor,
    Pade
};

namespace adaai {

    template<typename F, MethodE M = MethodE::Pade>
    constexpr F exp(F x) {  // return e^x
        static_assert(std::is_floating_point<F>::value, "Not a floating point number!");

        /*
            We can put long double instead of F in the intermediate
            calculations to achieve maximum accuracy, but in this case
            we may lose the speed of calculations
        */
        if (std::isnan(x)) {
            return x;
        }
        F y = x / C_LN_2<F>;  // e^x = 2^(x/ln2) = 2^y
        F y_int, y_float;
        y_float = std::modf(y, &y_int);

        if (y_int < INT32_MIN) {
            return 0.0;
        }
        if (y_int > INT32_MAX) {
            return INFINITY;
        }
        if (std::abs(y_float) > 0.5) {
            y_int += (std::signbit(y_float) ? -1 : 1);
            y_float += (std::signbit(y_float) ? 1 : -1);
        }
        F x1 = y_float * C_LN_2<F>; // e^x = 2^[x] * e^{{x}*Ln2}. x1 = {x}*Ln2
        F f1 = 0.0; // want find e^x1
        if constexpr (M == MethodE::Taylor) {
            F next = 1.0;
            unsigned long long n = 0;
            unsigned precl_n = exp_consts::MKExpTaylorOrder<F>();
            for (unsigned i = 0; i < precl_n; ++i) {
                f1 += next;
                next *= x1;
                next /= ++n;
            }
        } else if constexpr (M == MethodE::Pade) {
//            std::vector<F> members = {30240, 15120, 3360, 420, 30, 1}; // for x^0, x^1, ..., x^5.
//            std::vector<F> members = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0,
//                                      90.0, 1.0}; // to x^9
            std::vector<F> members = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1}; // to x^7
            F x_st = 1;
            F divisible = 0; // делимое
            F divisor = 0; // делитель
            for (unsigned st = 0; st < members.size(); ++st) {
                F number = x_st * members[st];
                divisible += number;
                if (st & 1) {
                    divisor -= number;
                } else {
                    divisor += number;
                }
                x_st *= x1;
            }
            f1 = divisible / divisor;
        }

        return std::ldexp(f1, int(y_int)); // f1 * 2^n
    }

} // namespace adaai

#endif // EXP_HPP_
