#ifndef EXP_TAYLOR_HPP_
#define EXP_TAYLOR_HPP_

#include "consts.hpp"

namespace adaai {

    namespace {
        template<typename F>
        constexpr inline unsigned MKExpTaylorOrder() {
            F ln2_n_sqrt2 = C_SQRT2<F>;
            F n_fact = 1;
            for (unsigned n = 0; n < 1000; n++) {
                if (ln2_n_sqrt2 < C_EPS<F> * n_fact) {
                    return n - 1;
                }
                ln2_n_sqrt2 *= C_LN2<F>;
                n_fact *= n + 1;
            }
            return 1000; // static_assert(false)
        }
    }

    template<typename F>
    constexpr F exp_taylor(F x) {
        F f = 0.0;
        F next = 1.0;
        unsigned n = 0;

        for (unsigned i = 0; i < MKExpTaylorOrder<F>(); ++i) {
            f += next;
            next *= x;
            next /= ++n;
        }

        return f;
    }

} // namespace adaai

#endif // EXP_TAYLOR_HPP_
