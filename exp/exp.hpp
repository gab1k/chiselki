#ifndef EXP_HPP_
#define EXP_HPP_

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "consts.hpp"
#include "exp_taylor.hpp"
#include "exp_pade.hpp"
#include "exp_chebyshev.hpp"

namespace adaai {

    enum class MethodE {
        Taylor,
        Pade,
        Chebyshev,
    };

    template<typename F, MethodE M = MethodE::Taylor>
    constexpr F exp(F x) {
        static_assert(std::is_floating_point_v<F>, "Not a floating point number");
        if (std::isnan(x)) {
            return x;
        }

        /*
         * e^x = 2^(log2(e) * x) = 2^(x/ln2) =
         * = 2^y = 2^(y_int + y_frac) =
         * = 2^y_int * e^(y_frac * ln2) =
         * = 2^y_int * e^x1
         * |x1| < 0.5 => Calculate for f1 = e^x1
         */

        F y = x / C_LN2<F>;
        F y_int, y_frac;
        y_frac = std::modf(y, &y_int);

        if (y_int < INT32_MIN) {
            return 0.0;
        }
        if (y_int > INT32_MAX) {
            return INFINITY;
        }

        if (std::abs(y_frac) > 0.5) {
            y_int += (std::signbit(y_frac) ? -1 : 1);
            y_frac += (std::signbit(y_frac) ? 1 : -1);
        }

        F x1 = y_frac * C_LN2<F>;

        F f1 = NAN;
        if constexpr (M == MethodE::Taylor) {
            f1 = exp_taylor(x1);
        } else if constexpr (M == MethodE::Pade) {
            f1 = exp_pade(x1);
        } else if constexpr (M == MethodE::Chebyshev) {
            f1 = exp_chebyshev(x1);
        } else {
            throw std::invalid_argument("Unexpected method");
        }

        return std::ldexp(f1, static_cast<uint32_t>(y_int));
    }

} // namespace adaai

#endif // EXP_HPP_
