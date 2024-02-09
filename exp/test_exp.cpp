#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "consts.hpp"
#include "exp.hpp"
#include <cmath>

template<typename F>
void checkExp(F value, F eps) {
    const F current = adaai::exp<MethodE::Pade, F>(value);
    const F expected = std::exp(value);
//    std::cout << current - expected << "\n";
    CHECK((std::abs(current - expected) <= eps ||
           std::isinf(expected) && std::isinf(current) ||
           std::isnan(expected) && std::isnan(current)));
}

template<typename F>
void testBasicWithTemplate() {
    F eps = adaai::C_EPS<F>;

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
