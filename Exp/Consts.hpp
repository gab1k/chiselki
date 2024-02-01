#include <cmath>
#include <cfloat>

namespace ADAAI {
// --------------------------------------- Define Ln2
    template<typename F>
    constexpr inline F
    Ln2;
    template<>
    inline float Ln2<float> = logf(2.0f);
    template<>
    inline double Ln2<double> = log(2.0);
    template<>
    inline long double Ln2<long double> = logl(2.0l);
// ---------------------------------------


// --------------------------------------- Define Eps
    template<typename F>
    constexpr inline F
    Eps;
    template<>
    constexpr inline float Eps<float> = FLT_EPSILON;
    template<>
    constexpr inline double Eps<double> = DBL_EPSILON;
    template<>
    constexpr inline long double Eps<long double> = LDBL_EPSILON;
// ---------------------------------------

// --------------------------------------- Define Sqrt2
    template<typename F>
    constexpr inline F
    Sqrt2;
    template<>
    constexpr inline float Sqrt2<float> = M_SQRT1_2f;
    template<>
    constexpr inline double Sqrt2<double> = M_SQRT1_2;
    template<>
    constexpr inline long double Sqrt2<long double> = M_SQRT1_2l;
// ---------------------------------------

}