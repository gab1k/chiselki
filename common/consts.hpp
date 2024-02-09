#ifndef CONST_HPP_
#define CONST_HPP_

#include <cfloat>
#include <cmath>

namespace adaai {
    // --------------------------------------- define C_LN_2
    template<typename F>
    constexpr inline F C_LN_2; // should be constexpr

    template<>
    constexpr inline float C_LN_2<float> = 1 / M_LOG2Ef;

    template<>
    constexpr inline double C_LN_2<double> = 1 / M_LOG2E;

    template<>
    constexpr inline long double C_LN_2<long double> = 1 / M_LOG2El;
    // --------------------------------------- define C_LN_2

    // --------------------------------------- define C_EPS
    template<typename F>
    constexpr inline F C_EPS;

    template<>
    constexpr inline float C_EPS<float> = FLT_EPSILON;

    template<>
    constexpr inline double C_EPS<double> = DBL_EPSILON;

    template<>
    constexpr inline long double C_EPS<long double> = LDBL_EPSILON;
    // --------------------------------------- define C_EPS


    // --------------------------------------- define C_SQRT_2
    template<typename F>
    constexpr inline F C_SQRT_2;

    template<>
    constexpr inline float C_SQRT_2<float> = M_SQRT1_2f;

    template<>
    constexpr inline double C_SQRT_2<double> = M_SQRT1_2;

    template<>
    constexpr inline long double C_SQRT_2<long double> = M_SQRT1_2l;
    // --------------------------------------- define C_SQRT_2

    namespace exp_consts {

        template<typename F>
        constexpr inline unsigned MKExpTaylorOrder() {
            F ln2_n_sqrt2 = C_SQRT_2<F>;
            F n_fact = 1;
            for (unsigned n = 0; n < 1000; n++) {
                if (ln2_n_sqrt2 / n_fact < C_EPS<F>) {
                    return n - 1;
                }
                ln2_n_sqrt2 *= C_LN_2<F>;
                n_fact *= n + 1;
            }
            return 1000; // static_assert(false)
        }

//        template<typename F>
//        constexpr inline std::vector<F> members_Pade{
//
//        };
    }
} // namespace adaai

#endif // CONST_HPP_