#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "consts.hpp"
#include "doctest.h"
#include "exp.hpp"

#include <random>
#include <cmath>


template<typename F, MethodE M = MethodE::Taylor>
void checkExp(F value, F eps) {
    const F expected = std::exp(value);
//    std::cout << get_eps<F>(expected) << "\n";
    const F current = adaai::exp<F, M>(value);
    CHECK((std::abs(current - expected) <= eps ||
           std::isinf(expected) && std::isinf(current) ||
           std::isnan(expected) && std::isnan(current)));
}

template<typename F, MethodE M = MethodE::Taylor>
void stress_test_exp(unsigned n, F from, F to) {
    srand(time(nullptr));
    std::cout.precision(100);
    for (unsigned i = 0; i < n; ++i) {
        F x = (F) rand() / RAND_MAX * (to - from) + from;
        checkExp<F, M>(x, adaai::C_EPS<F> * 4);
    }
}

template<typename F>
void testBasicWithTemplate() {
    F eps = adaai::C_EPS<F>;
    stress_test_exp<F, MethodE::Taylor>(100000, -0.34, 0.34);
    SUBCASE("Zero") {
        checkExp<F>(0, eps);
    }

    SUBCASE("Inf") {
        checkExp<F>(12345.1234, eps);
        checkExp<F>(INFINITY, eps);
        checkExp<F>(-INFINITY, eps);
    }

    SUBCASE("NaN") {
        checkExp<F>(NAN, eps);
    }

    SUBCASE("Small") {
        checkExp<F>(0.1234321, eps);
        checkExp<F>(0.0043212, eps);
        checkExp<F>(0.0000392, eps);
        checkExp<F>(0.0000087, eps);
        checkExp<F>(0.0000007, eps);
    }

    SUBCASE("Large") {
        eps = 1;
        checkExp<F>(1.123423, eps);
        checkExp<F>(2.432132, eps);
        checkExp<F>(10.09876, eps);
        checkExp<F>(13.33212, eps);
        checkExp<F>(15.49921, eps);
    }
}


TEST_CASE("Float") {
    testBasicWithTemplate<float>();
}

TEST_CASE("Double") {
    testBasicWithTemplate<double>();

    SUBCASE("Small extra case") {
        double eps = adaai::C_EPS<double>;
        checkExp<double>(0.0000000987645657, eps);
        checkExp<double>(0.0000000003213696, eps);
        checkExp<double>(0.0000000000013326, eps);
        checkExp<double>(0.0000000000000696, eps);
        checkExp<double>(0.0000000000000006, eps);
    }

    SUBCASE("Large extra case") {
        double eps = 1.;
        checkExp<double>(21.58393, eps);
        checkExp<double>(28.43232, eps);
        checkExp<double>(34.03276, eps);
    }
}

TEST_CASE("Long double") {
    testBasicWithTemplate<long double>();

    SUBCASE("Small extra case") {
        long double eps = adaai::C_EPS<long double>;
        checkExp<long double>(0.0000000000987631579, eps);
        checkExp<long double>(0.0000000000000000672, eps);
        checkExp<long double>(0.0000000000000000004, eps);
    }

    SUBCASE("Large extra case") {
        long double eps = 1.;
        checkExp<long double>(35.9274, eps);
        checkExp<long double>(39.423482, eps);
        checkExp<long double>(42.032026, eps);
    }
}

//TEST_CASE("Check func") {
//    stress_test_exp<long double>(10, 2, 10);
//}