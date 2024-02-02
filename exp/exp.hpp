#ifndef EXP_HPP_
#define EXP_HPP_

#include "consts.hpp"
#include <cstdint>
#include <cmath>

namespace adaai {

    template<typename F>
    constexpr F exp(F x) {  // return e^x
        static_assert(std::is_floating_point<F>::value, "Not a floating point number!");

        F y = x / C_LN_2<F>;  // e^x = 2^(x/ln2) = 2^y
        F y_int, y_float;
        y_float = std::modf(y, &y_int);

        if (y_int < INT32_MIN) {
            return 0.0;
        } else if (y_int > INT32_MAX) {
            return INFINITY;
        }

        if(y_float > 0.5){  // |y_float| should be < 0.5
            y_int ++;
            y_float -= 1.0;
        }

        if(y_float < -0.5){  // |y_float| should be < 0.5
            y_int --;
            y_float += 1.0;
        }

        F x1 = y_float * C_LN_2<F>; // e^x = 2^[x] * e^{{x}*Ln2}. x1 = {x}*Ln2
        F eps = C_EPS<F> * 10.0; // error rate

        // Start Taylor with n = 1, because 1st summand = 1
        F f1 = 1; // f1 = e^x1 - count by Taylor
        int n = 1;
        F x_n = x1; // x^n
        F fact_n = 1; // n!

        while (std::abs(C_SQRT_2<F> * x_n / fact_n) >= eps) {
            f1 += x_n / fact_n;
            n++;
            x_n *= x1;
            fact_n *= n;
        }

        return std::ldexp(f1, int(y_int)); // f1 * 2^n
    }

} // namespace adaai

#endif // EXP_HPP_
