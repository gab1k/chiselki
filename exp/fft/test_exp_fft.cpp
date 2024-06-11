#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "fft/exp_fft.hpp"

namespace adaai {

    TEST_CASE("FFT") {
        std::string filePath = "../exp/fft/logs/fft_logs.txt";
        std::ofstream ofs(filePath.c_str(), std::ios_base::out);
        ofs << std::setprecision(10) << "Testing float:\n";
        exp_fft<float>(ofs);
        ofs << "\n\nTesting double:\n";
        exp_fft<double>(ofs);
        ofs << "\n\nTesting long double:\n";
        exp_fft<long double>(ofs);
        ofs.close();
    }

} // namespace adaai
