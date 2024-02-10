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
    F y = num / adaai::C_LN_2<F>;
    F y_int;
    std::modf(y, &y_int);
    y_int = abs(int(y_int)) + 1;
    // т.к. в процессе решения мы можем поменять y_int, то максимум увеличим на 1. Это и будет максимальная ошибка
    if (M == MethodE::Taylor) {
        return std::ldexp(adaai::C_EPS<F> * 5, y_int); // C_EPS<F> * 5 - ошибка при |x| < 0.34 (ln2/2)
    } else if (M == MethodE::Pade) {
        F ln_2_st = 1;
        for (unsigned n = 1; n <= 19; ++n) { // n должен ходить до n + m + 1, где n,m - степени членов многочленов паде
            ln_2_st *= adaai::C_LN_2<F>;
        }
        return std::ldexp(ln_2_st, y_int); // ошибка O(x^{n+m+1})
    }
    return 0;
}

template<typename F, MethodE M = MethodE::Pade>
std::pair<F, F> checkExp(F value) {
    const F expected = std::exp(value);
    const F current = adaai::exp<F, M>(value);
    F eps = get_eps_exp<F, M>(value);
    CHECK((std::abs(current - expected) <= eps ||
           (std::isinf(expected) && std::isinf(current) && (std::signbit(expected) == std::signbit(current))) ||
           std::isnan(expected) && std::isnan(current)));
    F abs_err = fabsl(expected - current);
    F rel_err = fabsl(current / expected - 1.0);
    return {abs_err, rel_err};
}

template<typename F, MethodE M = MethodE::Taylor>
std::pair<F, F> stress_test_exp(unsigned n, F from, F to) {
    srand(time(nullptr));
    std::pair<F, F> max_errs = {0, 0};
    for (unsigned i = 0; i < n; ++i) {
        F x = (F) rand() / RAND_MAX * (to - from) + from;
        std::pair<F, F> errs = checkExp<F, M>(x);
        max_errs = {std::max(errs.first, max_errs.first), std::max(errs.second, max_errs.second)};
    }
    return max_errs;
}

template<typename F>
void testBasicWithTemplate() {

    SUBCASE("Close to ln2/2") {
        checkExp<F>(-0.3465735902);
        checkExp<F>(0.3465735902);
    }

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

template<typename F, MethodE M = MethodE::Taylor>
void log_info(std::ofstream &ofs, F min_value, F max_value) {
    std::pair<F, F> errs_p = stress_test_exp<F, M>(10000, min_value, max_value);
    std::pair<F, F> errs_m = stress_test_exp<F, M>(10000, -max_value, -min_value);

    std::pair<F, F> max_errs = {std::max(errs_m.first, errs_p.first), std::max(errs_m.second, errs_p.second)};
    ofs << "Case: " << min_value << " < |x| < " << max_value << "\n";
    ofs << "Absolute error = " << max_errs.first << "; It is " << max_errs.first / adaai::C_EPS<F> << " * eps\n";
    ofs << "Relative error = " << max_errs.second << "; It is " << max_errs.second / adaai::C_EPS<F> << " * eps\n";
    F our_err = std::max(errs_p.second, errs_m.first);
    ofs << "Our error is = " << our_err << "; It is " << our_err / adaai::C_EPS<F> << " * eps\n\n";

}

TEST_CASE("Log Error Float") {
    std::string filePath = "../exp/float_logs.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out);
    ofs << "Eps = " << adaai::C_EPS<float> << "\n\n";
    ofs << "Taylor:\n";
    float hod[] = {0, 0.34, 3.0, 5.0, 7.0};
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<float, MethodE::Taylor>(ofs, hod[i - 1], hod[i]);
    }

    ofs << "\n\nPade:\n";
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<float, MethodE::Pade>(ofs, hod[i - 1], hod[i]);
    }

    ofs.close();
}

TEST_CASE("Log Error Double") {
    std::string filePath = "../exp/double_logs.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out);
    ofs << "Eps = " << adaai::C_EPS<double> << "\n\n";
    ofs << "Taylor:\n";

    double hod[] = {0.00, 0.34, 3.00, 5.00, 7.00, 10.0, 15.0};
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<double, MethodE::Taylor>(ofs, hod[i - 1], hod[i]);
    }

    ofs << "\n\nPade:\n";
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<double, MethodE::Pade>(ofs, hod[i - 1], hod[i]);
    }

    ofs.close();
}

TEST_CASE("Log Error Long Double") {
    std::string filePath = "../exp/long_double_logs.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out);
    ofs << "Eps = " << adaai::C_EPS<long double> << "\n\n";
    ofs << "Taylor:\n";
    long double hod[] = {0, 0.34, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0};
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<long double, MethodE::Taylor>(ofs, hod[i - 1], hod[i]);
    }

    ofs << "\n\nPade:\n";
    for (int i = 1; i < sizeof(hod) / sizeof(hod[0]); ++i) {
        log_info<long double, MethodE::Pade>(ofs, hod[i - 1], hod[i]);
    }
    ofs.close();
}