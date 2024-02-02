#ifndef EXP_HPP_
#define EXP_HPP_

#include "consts.hpp"
#include <cmath>
#include <cstdint>

namespace adaai {

    template<typename F>
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

        F f1 = 0.0; // f1 = e^x1 - count by Taylor
        F next = 1.0;
        unsigned long long n = 0;

        while (std::abs(next * C_SQRT_2<F>) >= C_EPS<F>) {
            f1 += next;
            next *= x1;
            next /= ++n;
        }

        return std::ldexp(f1, int(y_int)); // f1 * 2^n
    }

} // namespace adaai

#endif // EXP_HPP_
