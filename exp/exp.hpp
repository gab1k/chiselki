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

    template<MethodE M = MethodE::Pade, typename F>
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
            unsigned precl_n = MKExpTaylorOrder<F>();
            for(unsigned i = 0; i < precl_n; ++i) {
                f1 += next;
                next *= x1;
                next /= ++n;
            }
//            while (std::abs(next * C_SQRT_2<F>) >= C_EPS<F>) {
//                f1 += next;
//                next *= x1;
//                next /= ++n;
//            }

            std::cout << n << "\n";
        } else if constexpr (M == MethodE::Pade) {
            std::vector<F> members = {30240, 15120, 3360, 420, 30, 1}; // for x^0, x^1, ..., x^4.
            F x_st = 1;
            F divisible = 0; // делимое
            F divisor = 0; // делитель
            for (unsigned st = 0; st <= 6; ++st) {
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
