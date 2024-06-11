#ifndef EXP_FFT_HPP_
#define EXP_FFT_HPP_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "consts.hpp"

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

namespace adaai {

    template<typename F>
    constexpr F exp_fft(std::ofstream &ofs) {
        int N = 32;
        ofs << "Start testing with N = " << 32 << "\n\n";
        double a[] = {14.089543728627603 / 2, 0, -7.678045762940514, 0, 2.8133624539524732, 0, -1.5306270816125131, 0,
                      0.82339375927268074, 0, -0.58479237046751483, 0, 0.37511808822149146, 0, -0.30093848048476058, 0,
                      0.21087657800537977, 0, -0.18081376531146307, 0, 0.13339960844918627, 0, -0.11913778714815918, 0,
                      0.090818445485041704, 0, -0.083286819158978232, 0, 0.064872286930420463, 0, -0.060540956229970837,
                      0, 0.047820151096467916, 0, -0.045117560882689439, 0, 0.035927269512624666, 0,
                      -0.034079659693320645, 0, 0.027209431900052818, 0, -0.025803163212257383, 0, 0.020529010894549, 0,
                      -0.019325460880796904, 0, 0.015190400509132839, 0, -0.014040122413933981, 0, 0.010742311490565823,
                      0, -0.0095418335808245976, 0, 0.006873627011116772, 0, -0.0055417595296704336, 0,
                      0.0033544710557849812, 0, -0.001817772246759458, 0};

        gsl_fft_complex_radix2_forward(a, 1, N);
        std::cout.precision(17);
        for (int i = 0; i < N; i += 2) {
            F xi = (Pi<F> * i) / N;
            ofs << "for i = " << i << " FFT value is: " << REAL(a, i / 2) << " expected: " << exp(xi) << "\n";
            ofs << "absolute error is: " << abs(exp(xi) - REAL(a, i / 2)) << "\n\n";
        }
        return 0;
    }

} // namespace adaai

#endif // EXP_FFT_HPP_
