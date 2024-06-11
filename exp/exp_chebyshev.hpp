#ifndef EXP_CHEBYCHEV_HPP_
#define EXP_CHEBYCHEV_HPP_

#include <gsl/gsl_chebyshev.h>

namespace adaai {

    template<typename F>
    constexpr F exp_chebyshev(F x) {
        int N = 20;
        auto *c = (double *) calloc(N + 1, sizeof(double));
        double d[] = {1.2660658777520084 * 2, 1.1303182079849703, 0.27149533953407662, 0.044336849848663783,
                      0.0054742404420937332, 0.00054292631191394389, 4.4977322954295149e-05, 3.1984364624019909e-06,
                      1.992124806672796e-07, 1.1036771725517344e-08, 5.5058960796737464e-10, 2.4979566169849831e-11,
                      1.0391522306785704e-12, 3.9912633564144015e-14, 1.4237580108256574e-15, 4.7409261025614987e-17,
                      1.4801800572077666e-18, 4.3499194966456986e-20, 1.2074283482289013e-21, 3.1774430216550031e-23,
                      7.9436075541375073e-25};
        for (int i = 0; i < N + 1; i++) {
            c[i] = d[i];
        }

        gsl_cheb_series *cs = gsl_cheb_alloc(N);
        cs->c = c;
        cs->a = -1;
        cs->b = 1;

        F f = gsl_cheb_eval(cs, x);
        gsl_cheb_free(cs);

        return f;
    }

} // namespace adaai

#endif // EXP_CHEBYCHEV_HPP_
