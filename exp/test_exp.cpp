#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "consts.hpp"
#include "doctest.h"
#include "exp.hpp"

#include <random>
#include <cmath>

template<typename F, MethodE M = MethodE::Taylor>
F get_eps_exp(F num) {
    if (std::isnan(num) || std::isinf(num)) {
        return 0;
    }
    if (M == MethodE::Taylor) {
        F y = num / adaai::C_LN_2<F>;
        F y_int;
        std::modf(y, &y_int);
        return std::ldexp(adaai::C_EPS<F> * 5, abs(int(y_int)) + 1); // C_EPS<F> * 5 - ошибка при |x| < 0.34 (ln2/2)
        // т.к. в процессе решения мы можем поменять y_int, то максимум увеличим на 1. Это и будет максимальная ошибка
    } else if (M == MethodE::Pade) {
        F ln_2_st = 1;
        for (unsigned n = 1; n <= 15; ++n) { // n должен ходить до n + m + 1, где n,m - степени членов многочленов паде
            ln_2_st *= adaai::C_LN_2<F>;
        }
        return ln_2_st;
    }
    return 0;
}

template<typename F, MethodE M = MethodE::Taylor>
void checkExp(F value) {
    const F expected = std::exp(value);
    const F current = adaai::exp<F, M>(value);
    F eps = get_eps_exp<F>(value);
    CHECK((std::abs(current - expected) <= eps ||
           (std::isinf(expected) && std::isinf(current) && (std::signbit(expected) == std::signbit(current))) ||
           std::isnan(expected) && std::isnan(current)));
}

template<typename F, MethodE M = MethodE::Taylor>
void stress_test_exp(unsigned n, F from, F to) {
    srand(time(nullptr));
    for (unsigned i = 0; i < n; ++i) {
        F x = (F) rand() / RAND_MAX * (to - from) + from;
        checkExp<F, M>(x);
    }
}

template<typename F>
void testBasicWithTemplate() {
    stress_test_exp<F, MethodE::Taylor>(100000, -5.34, 5.34);
    SUBCASE("Zero") {
        checkExp<F>(0);
    }

    SUBCASE("Inf") {
        checkExp<F>(12345.1234);
        checkExp<F>(INFINITY);
        checkExp<F>(-INFINITY);
    }

    SUBCASE("NaN") {
        checkExp<F>(NAN);
    }

    SUBCASE("Small") {
        checkExp<F>(0.1234321);
        checkExp<F>(0.0043212);
        checkExp<F>(0.0000392);
        checkExp<F>(0.0000087);
        checkExp<F>(0.0000007);
    }

    SUBCASE("Large") {
        checkExp<F>(1.123423);
        checkExp<F>(2.432132);
        checkExp<F>(10.09876);
        checkExp<F>(13.33212);
        checkExp<F>(15.49921);
    }
}


TEST_CASE("Float") {
    testBasicWithTemplate<float>();
}

TEST_CASE("Double") {
    testBasicWithTemplate<double>();

    SUBCASE("Small extra case") {
        checkExp<double>(0.0000000987645657);
        checkExp<double>(0.0000000003213696);
        checkExp<double>(0.0000000000013326);
        checkExp<double>(0.0000000000000696);
        checkExp<double>(0.0000000000000006);
    }

    SUBCASE("Large extra case") {
        checkExp<double>(21.58393);
        checkExp<double>(28.43232);
        checkExp<double>(34.03276);
    }
}

TEST_CASE("Long double") {
    testBasicWithTemplate<long double>();

    SUBCASE("Small extra case") {
        checkExp<long double>(0.0000000000000031579);
        checkExp<long double>(0.0000000000000000672);
        checkExp<long double>(0.0000000000000000004);
    }

    SUBCASE("Large extra case") {
        checkExp<long double>(35.9274);
        checkExp<long double>(39.423482);
        checkExp<long double>(42.032026);
    }
}

//TEST_CASE("Check func") {
//    stress_test_exp<long double>(10, 2, 10);
//}