#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "consts.hpp"
#include "exp.hpp"
#include <iostream>

template<typename F>
void checkExp(F value, F eps, bool debug=false) {
    const F current = adaai::exp<F>(value);
    const F expected = std::exp(value);
    CHECK((std::abs(current - expected) < eps ||
            std::isinf(expected) && std::isinf(current) ||
            std::isnan(expected) && std::isnan(current)));
    if (debug) {
        std::cout.precision(20);
        std::cout << "Value: " << value << "\n" <<
                    "Expected: " << expected << "\n" <<
                    "Current: " << current << "\n";
    }
}

template<typename F>
void testWithTemplate() {
    SUBCASE("Small") {
        checkExp<F>(0., 0.00001);
        checkExp<F>(0.0005, 0.00001);
        checkExp<F>(0.123456789, 0.00001);
        checkExp<F>(0.00000000005, 0.01);
        checkExp<F>(0.00000000000000005, 0.01);
    }
    SUBCASE("Large") {
        checkExp<F>(12., 0.01, true);
    }
    SUBCASE("Inf") {
        checkExp<F>(105.5, 0.01);
        checkExp<F>(112305.5, 0.01);
    }
    SUBCASE("Nan") {
        checkExp<F>(NAN, 0.01);
    }
}

TEST_CASE("Test exp function") {
    SUBCASE("Float") {
        testWithTemplate<float>();
    }
    SUBCASE("Double") {
        testWithTemplate<float>();
    }
    SUBCASE("Long double") {
        testWithTemplate<float>();
    }
}
