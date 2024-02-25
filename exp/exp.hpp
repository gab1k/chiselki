#ifndef EXP_HPP_
#define EXP_HPP_

#include "consts.hpp"
#include <cmath>
#include <cstdint>

#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

enum class MethodE : int {
    Taylor,
    Pade,
    Chebyshev
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
            constexpr unsigned precl_n = exp_consts::MKExpTaylorOrder<F>();
            for (unsigned i = 0; i < precl_n; ++i) {
                f1 += next;
                next *= x1;
                next /= ++n;
            }
        } else if constexpr (M == MethodE::Pade) {
//            std::vector<F> members = {30240, 15120, 3360, 420, 30, 1}; // for x^0, x^1, ..., x^5.
//            std::vector<F> members = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1}; // to x^7
            std::vector<F> members = {17643225600.0, 8821612800.0, 2075673600.0, 302702400.0, 30270240.0, 2162160.0, 110880.0, 3960.0,
                                      90.0, 1.0}; // to x^9
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
        } else if constexpr (M == MethodE::Chebyshev) {
            int N = 14;
            auto *c = (double *)calloc(N + 1, sizeof(double));
            double d[] = {1.61186949110557,
                          1.2996917014510012,
                          0.62435557930913788,
                          0.10196108566545124,
                          0.012589065316430418,
                          0.0012485631340078945,
                          0.00010343397635147418,
                          7.3554177902041539e-06,
                          4.5812728861601064e-07,
                          2.5381172347981936e-08,
                          1.2661863523359072e-09,
                          5.7445301263787905e-11,
                          2.3897245325735763e-12,
                          9.1912482022060619e-14,
                          3.2825886436450219e-15};
            for(int i = 0; i < N + 1; i ++){
                c[i] = d[i];
            }
            gsl_cheb_series *cs = gsl_cheb_alloc(N);
            cs->c = c;
            f1 = gsl_cheb_eval(cs, x1);
            gsl_cheb_free(cs);
        }

        return std::ldexp(f1, int(y_int)); // f1 * 2^n
    }

} // namespace adaai

#endif // EXP_HPP_
