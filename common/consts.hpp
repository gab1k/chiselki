#ifndef CONSTS_HPP_
#define CONSTS_HPP_

#include <cfloat>
#include <cmath>

namespace adaai {

    template<typename F>
    constexpr inline F C_LN2;
    template<>
    constexpr inline float C_LN2<float> = 1 / M_LOG2Ef;
    template<>
    constexpr inline double C_LN2<double> = 1 / M_LOG2E;
    template<>
    constexpr inline long double C_LN2<long double> = 1 / M_LOG2El;

    template<typename F>
    constexpr inline F C_EPS;
    template<>
    constexpr inline float C_EPS<float> = FLT_EPSILON;
    template<>
    constexpr inline double C_EPS<double> = DBL_EPSILON;
    template<>
    constexpr inline long double C_EPS<long double> = LDBL_EPSILON;

    template<typename F>
    constexpr inline F C_SQRT2;
    template<>
    constexpr inline float C_SQRT2<float> = M_SQRT1_2f;
    template<>
    constexpr inline double C_SQRT2<double> = M_SQRT1_2;
    template<>
    constexpr inline long double C_SQRT2<long double> = M_SQRT1_2l;

    namespace exp_consts {

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

    // --------------------------------------- define Pi
    template<typename F>
    constexpr inline F Pi;

    template<>
    constexpr inline float Pi<float> = 3.1415926535f;

    template<>
    constexpr inline double Pi<double> = 3.14159265358979323846;

    template<>
    constexpr inline long double Pi<long double> = 3.1415926535897932384626433846l;
    // --------------------------------------- define Pi
} // namespace adaai

#endif // CONSTS_HPP_