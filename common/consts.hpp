#ifndef CONST_HPP_
#define CONST_HPP_

#include <cfloat>
#include <cmath>

namespace adaai {
    // --------------------------------------- define C_LN_2
    template<typename F>
    constexpr inline F C_LN_2;

    template<>
    inline float C_LN_2<float> = logf(2.0f);

    template<>
    inline double C_LN_2<double> = log(2.0);

    template<>
    inline long double C_LN_2<long double> = logl(2.0l);
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
} // namespace adaai

#endif // CONST_HPP_