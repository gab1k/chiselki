#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "consts.hpp"
#include "exp.hpp"

template<typename F>
F absDifferent(F value) {
    return std::abs(adaai::exp(value) - std::exp(value));
}

template<typename F>
void testWithTemplate() {
    CHECK((absDifferent<F>(0.) < adaai::C_EPS<F>));
    CHECK((absDifferent<F>(0.5) < adaai::C_EPS<F>));
    CHECK((absDifferent<F>(1.5) < adaai::C_EPS<F>));
}

TEST_CASE("Test exp function") {
    SUBCASE("Test float") {
        testWithTemplate<float>();
    }
    SUBCASE("Test double") {
        testWithTemplate<float>();
    }
    SUBCASE("Test long double") {
        testWithTemplate<float>();
    }
}
